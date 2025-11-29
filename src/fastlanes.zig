//! This file implements bit packing and delta coding.
//!
//! Code is a port of fastlanes Rust implementation by spiraldb:
//! github.com/spiraldb/fastlanes
//!
//! It was specifically ported from commit hash 8b655cf.
//!
//! Same license as the original code can be found in project root `LICENSE-APACHE`

const std = @import("std");

const FL_ORDER = [_]usize{ 0, 4, 2, 6, 1, 5, 3, 7 };

fn transpose_idx(i: usize) usize {
    const lane = i % 16;
    const order = (i / 16) % 8;
    const row = i / 128;

    return (lane * 64) + (FL_ORDER[order] * 8) + row;
}

fn index(bit_idx: usize, lane: usize) usize {
    const o = bit_idx / 8;
    const s = bit_idx % 8;
    return (FL_ORDER[o] * 16) + (s * 128) + lane;
}

pub fn FastLanes(comptime T: type) type {
    return struct {
        pub const N_BITS = @sizeOf(T) * 8;
        pub const N_LANES = 1024 / N_BITS;

        fn mask(comptime width: comptime_int) T {
            return comptime (1 << width) - 1;
        }

        fn transpose(input: [1024]T) [1024]T {
            var output: [1024]T = undefined;
            for (0..1024) |i| {
                output[i] = input[transpose_idx(i)];
            }

            return output;
        }

        fn untranspose(input: [1024]T) [1024]T {
            var output: [1024]T = undefined;
            for (0..1024) |i| {
                output[transpose_idx(i)] = input[i];
            }

            return output;
        }

        pub fn Packer(comptime W: comptime_int) type {
            if (W > N_BITS) {
                @compileError("invalid W");
            }
            return struct {
                pub const PACKED_LEN = 1024 * W / N_BITS;

                pub fn pack(input: [1024]T) [PACKED_LEN]T {
                    @setEvalBranchQuota(4096);

                    if (W == 0) {
                        return .{};
                    }
                    if (W > N_BITS) {
                        @compileError("invalid W");
                    }

                    const MASK = mask(W);
                    var out: [PACKED_LEN]T = undefined;

                    for (0..N_LANES) |lane_idx| {
                        if (W == N_BITS) {
                            inline for (0..N_BITS) |bit_idx| {
                                const idx = index(bit_idx, lane_idx);
                                out[N_LANES * bit_idx + lane_idx] = input[idx];
                            }
                        } else {
                            var tmp: T = 0;
                            inline for (0..N_BITS) |bit_idx| {
                                const idx = index(bit_idx, lane_idx);
                                const src = input[idx] & MASK;

                                if (bit_idx == 0) {
                                    tmp = src;
                                } else {
                                    tmp |= src << (bit_idx * W) % N_BITS;
                                }

                                const curr_word = bit_idx * W / N_BITS;
                                const next_word = (bit_idx + 1) * W / N_BITS;

                                if (next_word > curr_word) {
                                    out[N_LANES * curr_word + lane_idx] = tmp;
                                    const remaining_bits = ((bit_idx + 1) * W) % N_BITS;
                                    tmp = src >> W - remaining_bits;
                                }
                            }
                        }
                    }

                    return out;
                }

                pub fn unpack(input: [PACKED_LEN]T) [1024]T {
                    @setEvalBranchQuota(4096);

                    if (W > N_BITS) {
                        @compileError("invalid W");
                    }

                    if (W == 0) {
                        var out: [1024]T = undefined;
                        for (0..1024) |i| {
                            out[i] = 0;
                        }
                        return out;
                    }

                    var out: [1024]T = undefined;

                    for (0..N_LANES) |lane_idx| {
                        if (W == N_BITS) {
                            inline for (0..N_BITS) |bit_idx| {
                                const idx = index(bit_idx, lane_idx);
                                const src = input[N_LANES * bit_idx + lane_idx];
                                out[idx] = src;
                            }
                        } else {
                            var src = input[lane_idx];
                            var tmp: T = 0;

                            inline for (0..N_BITS) |bit_idx| {
                                const curr_word = bit_idx * W / N_BITS;
                                const next_word = (bit_idx + 1) * W / N_BITS;

                                const shift = (bit_idx * W) % N_BITS;

                                if (next_word > curr_word) {
                                    const remaining_bits = ((bit_idx + 1) * W) % N_BITS;
                                    const current_bits = W - remaining_bits;
                                    tmp = (src >> shift) & mask(current_bits);

                                    if (next_word < W) {
                                        src = input[N_LANES * next_word + lane_idx];
                                        tmp |= (src & mask(remaining_bits)) << current_bits;
                                    }
                                } else {
                                    tmp = (src >> shift) & mask(W);
                                }

                                const idx = index(bit_idx, lane_idx);
                                out[idx] = tmp;
                            }
                        }
                    }

                    return out;
                }
            };
        }

        pub fn delta(input: [1024]T, base: [N_LANES]T) [1024]T {
            var out: [1024]T = undefined;

            for (0..N_LANES) |lane_idx| {
                var prev = base[lane_idx];

                inline for (0..N_BITS) |bit_idx| {
                    const idx = index(bit_idx, lane_idx);
                    const next = input[idx];
                    out[idx] = next -% prev;
                    prev = next;
                }
            }

            return out;
        }

        pub fn undelta(input: [1024]T, base: [N_LANES]T) [1024]T {
            var out: [1024]T = undefined;

            for (0..N_LANES) |lane_idx| {
                var prev = base[lane_idx];

                inline for (0..N_BITS) |bit_idx| {
                    const idx = index(bit_idx, lane_idx);
                    const next = input[idx] +% prev;
                    out[idx] = next;
                    prev = next;
                }
            }

            return out;
        }
    };
}

// test "pack to" {
//     var data: [1024]u64 = undefined;
//     for (0..1024) |i| {
//         data[i] = i;
//     }

//     const FL = FastLanes(u64);

//     const packed_d = FL.pack(10, data);

//     const unpacked_d = FL.unpack(10, packed_d);

//     for (0..1024) |i| {
//         std.debug.assert(unpacked_d[i] == data[i]);
//     }
// }

test "delta pack to" {
    const base = -213123213;

    var data: [1024]i64 = undefined;
    for (0..1024) |i| {
        const v: i64 = @intCast(i);
        data[i] = base + v * 2;
    }

    const FL = FastLanes(i64);
    const Packer = FL.Packer(2);

    const transposed_d = FL.transpose(data);

    const bases = transposed_d[0..FL.N_LANES].*;

    const delta_d = FL.delta(transposed_d, bases);

    const packed_d = Packer.pack(delta_d);

    const unpacked_d = Packer.unpack(packed_d);

    const undelta_d = FL.undelta(unpacked_d, bases);

    const untransposed_d = FL.untranspose(undelta_d);

    for (0..1024) |i| {
        std.debug.assert(untransposed_d[i] == data[i]);
    }
}
