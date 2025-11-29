const std = @import("std");

pub const Error = error{Invalid};

const FL_ORDER = [_]usize{ 0, 4, 2, 6, 1, 5, 3, 7 };

fn transpose(comptime T: type, input: [1024]T) [1024]T {
    var output: [1024]T = undefined;
    for (0..1024) |i| {
        output[i] = input[transpose_idx(i)];
    }

    return output;
}

fn untranspose(comptime T: type, input: [1024]T) [1024]T {
    var output: [1024]T = undefined;
    for (0..1024) |i| {
        output[transpose_idx(i)] = input[i];
    }

    return output;
}

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

fn FastLanes(comptime T: type) type {
    return struct {
        pub const N_BITS = @sizeOf(T) * 8;
        pub const N_LANES = 1024 / N_BITS;

        pub fn packed_len(comptime W: comptime_int) comptime_int {
            return 1024 * W / N_BITS;
        }

        fn mask(comptime width: comptime_int) T {
            return comptime (1 << width) - 1;
        }

        pub fn pack(comptime W: comptime_int, input: [1024]T) [packed_len(W)]T {
            @setEvalBranchQuota(4096);

            if (W > N_BITS) {
                @compileError("invalid W");
            }

            if (W == 0) {
                return .{};
            }

            const MASK = mask(W);
            var out: [packed_len(W)]T = undefined;

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

        pub fn unpack(comptime W: comptime_int, input: [packed_len(W)]T) [1024]T {
            @setEvalBranchQuota(4096);

            if (W > N_BITS) {
                @compileError("invalid W");
            }

            if (W == 0) {
                return std.mem.zeroes([1024]T);
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

// fn Zint(comptime T: type) type {
//     const FL = FastLanes(T);

//     return struct {
//         pub fn pack(input: [1024]T, noalias out: []u8) usize {
//             std.debug.assert(out.len >= 1 + @sizeOf(T) * 1024);

//             var max = input[0];
//             for (0..1024) |idx| {
//                 max = @max(input[idx], max);
//             }

//             const num_bits: u8 = std.math.log2_int_ceil(T, max);

//             out[0] = num_bits;
//             inline for (0..FL.N_BITS) |nb| {
//                 if (nb == num_bits) {
//                     const packed_input = FL.pack(nb, input);
//                     @as([][FL.packed_len(nb)]T, @ptrCast(out[1..1+nb*1024/8]))[0] = packed_input;
//                     return 1 + nb * 1024 / 8;
//                 }
//             }

//             std.debug.assert(num_bits == FL.N_BITS);
//             @memcpy(out[1..1+@sizeOf(T) * 1024], input);

//             return 1 + @sizeOf(T) * 1024;
//         }

//         pub fn unpack(noalias input: []const u8) Error!struct { [1024]T, usize } {
//             if (input.len == 0) {
//                 return Error.Invalid;
//             }
//             const num_bits = input[0];

//             inline for (0..FL.N_BITS) |nb| {
//                 if (nb == num_bits) {
//                     const input_end = 1 + FL.packed_len(nb) * @sizeOf(T);
//                     if (input.len < input_end) {
//                         return Error.Invalid;
//                     }
//                     const packed_input: [FL.packed_len(nb)]T = @bitCast(input[1..input_end]);
//                     return .{ FL.unpack(nb, packed_input), input_end };
//                 }
//             }

//             std.debug.assert(num_bits == FL.N_BITS);

//             const input_end = 1 + 1024 * @sizeOf(T);
//             if (input.len < input_end) {
//                 return Error.Invalid;
//             }
//             return .{ @bitCast(input[1..input_end]), input_end };
//         }
//     };
// }

test "pack to" {
    var data: [1024]u64 = undefined;
    for (0..1024) |i| {
        data[i] = i;
    }

    const FL = FastLanes(u64);

    const packed_d = FL.pack(10, data);

    const unpacked_d = FL.unpack(10, packed_d);

    for (0..1024) |i| {
        std.debug.assert(unpacked_d[i] == data[i]);
    }
}

test "delta pack to" {
    const base = 213123213;

    var data: [1024]u64 = undefined;
    for (0..1024) |i| {
        data[i] = @intCast(i + base);
    }

    const FL = FastLanes(u64);

    const transposed_d = transpose(u64, data);

    const bases = transposed_d[0..FL.N_LANES].*;

    const delta_d = FL.delta(transposed_d, bases);

    const packed_d = FL.pack(1, delta_d);

    const unpacked_d = FL.unpack(1, packed_d);

    const undelta_d = FL.undelta(unpacked_d, bases);

    const untransposed_d = untranspose(u64, undelta_d);

    for (0..1024) |i| {
        std.debug.assert(untransposed_d[i] == data[i]);
    }
}

// test "zint pack" {
//     var data: [1024]u64 = undefined;
//     for (0..1024) |i| {
//         data[i] = i;
//     }

//     const Z = Zint(u64);

//     var packed_d: [1024 * @sizeOf(u64)]u8 = undefined;
//     const packed_len = Z.pack(data, &packed_d);

//     const unpacked_d, const consumed = Z.unpack(&packed_d);

//     std.debug.assert(consumed == packed_len);

//     for (0..1024) |i| {
//         std.debug.assert(unpacked_d[i] == data[i]);
//     }
// }
