//! This file implements bit packing and delta coding.
//!
//! Code is a port of fastlanes Rust implementation by spiraldb:
//! github.com/spiraldb/fastlanes
//!
//! It was specifically ported from commit hash 8b655cf.
//!
//! Same license as the original code can be found in project root `LICENSE-APACHE`

const FL_ORDER = [_]usize{ 0, 4, 2, 6, 1, 5, 3, 7 };

inline fn index(row: usize, lane: usize) usize {
    const o = row / 8;
    const s = row % 8;
    return (FL_ORDER[o] * 16) + (s * 128) + lane;
}

inline fn transpose_idx(idx: usize) usize {
    const lane = idx % 16;
    const order = (idx / 16) % 8;
    const row = idx / 128;

    return (lane * 64) + (FL_ORDER[order] * 8) + row;
}

        // const N_BITS = @sizeOf(T) * 8;
        // pub const N_LANES = 1024 / N_BITS;

pub fn n_lanes(comptime T: type) comptime_int {
    return 1024 / @bitSizeOf(T);
}

        inline fn mask(comptime T: type, width: comptime_int) T {
            return (1 << width) - 1;
        }

        pub fn rle_encode(
            comptime T: type,
            noalias input: *const [1024]T,
            noalias rle_vals: *[1024]T,
            noalias rle_idxs: *[1024]u16,
        ) usize {
            rle_vals[0] = input[0];
            rle_idxs[0] = 0;

            var rle_val_idx: u16 = 0;
            var prev = input[0];
            for (1..1024) |i| {
                const curr = input[i];
                if (curr != prev) {
                    rle_val_idx += 1;
                    rle_vals[rle_val_idx] = curr;
                    prev = curr;
                }
                rle_idxs[i] = rle_val_idx;
            }

            return rle_val_idx + 1;
        }

        pub fn rle_decode(
            comptime T: type,
            noalias rle_vals: []const T,
            noalias rle_idxs: *const [1024]u16,
            noalias output: *[1024]T,
        ) void {
            for (rle_idxs, output) |i, *out| {
                out.* = rle_vals[i];
            }
        }

        pub fn delta(
            comptime T: type,
            noalias input: *const [1024]T,
            noalias base: *const [n_lanes(T)]T,
            noalias output: *[1024]T,
        ) void {
            for (0..n_lanes(T)) |lane| {
                var prev = base[lane];
                inline for (0..@bitSizeOf(T)) |row| {
                    const idx = index(row, lane);
                    const next = input[idx];
                    output[idx] = next -% prev;
                    prev = next;
                }
            }
        }

        pub fn undelta(
            comptime T: type,
            noalias input: *const [1024]T,
            noalias base: *const [n_lanes(T)]T,
            noalias output: *[1024]T,
        ) void {
            for (0..n_lanes(T)) |lane| {
                var prev = base[lane];
                inline for (0..@bitSizeOf(T)) |row| {
                    const idx = index(row, lane);
                    const next = input[idx] +% prev;
                    output[idx] = next;
                    prev = next;
                }
            }
        }

        pub fn transpose(
            comptime T: type,
            noalias input: *const [1024]T,
            noalias output: *[1024]T,
        ) void {
            const UNROLL = 32;
            for (0..1024 / UNROLL) |i| {
                const idx = i * UNROLL;
                for (0..UNROLL) |j| {
                    output[idx + j] = input[transpose_idx(idx + j)];
                }
            }
        }

        pub fn untranspose(
            comptime T: type,
            noalias input: *const [1024]T,
            noalias output: *[1024]T,
        ) void {
            const UNROLL = 32;
            for (0..1024 / UNROLL) |i| {
                const idx = i * UNROLL;
                for (0..UNROLL) |j| {
                    output[transpose_idx(idx + j)] = input[idx + j];
                }
            }
        }

        pub fn dyn_bit_pack(
            comptime T: type,
            noalias input: *const [1024]T,
            noalias output: []align(1) T,
            width: usize,
        ) usize {
            inline for (0..@bitSizeOf(T) + 1) |W| {
                if (W == width) {
                    @call(.never_inline, bit_pack, .{ T, W, input, output[0..packed_len(T, W)] });
                    return packed_len(T, W);
                }
            }
            unreachable;
        }

        pub fn dyn_bit_unpack(
            comptime T: type,
            noalias input: []align(1) const T,
            noalias output: *[1024]T,
            width: usize,
        ) usize {
            inline for (0..@bitSizeOf(T) + 1) |W| {
                if (W == width) {
                    @call(.never_inline, bit_unpack, .{ T, W, input[0..packed_len(T, W)], output });
                    return packed_len(T, W);
                }
            }
            unreachable;
        }

        pub fn dyn_for_pack(
            comptime T: type,
            noalias input: *const [1024]T,
            reference: T,
            noalias output: []align(1) T,
            width: usize,
        ) usize {
            inline for (0..@bitSizeOf(T) + 1) |W| {
                if (W == width) {
                    @call(.never_inline, for_pack, .{ T, W, input, reference, output[0..packed_len(T, W)] });
                    return packed_len(T, W);
                }
            }
            unreachable;
        }

        pub fn dyn_for_unpack(
            comptime T: type,
            noalias input: []align(1) const T,
            reference: T,
            noalias output: *[1024]T,
            width: usize,
        ) usize {
            inline for (0..@bitSizeOf(T) + 1) |W| {
                if (W == width) {
                    @call(.never_inline, for_unpack, .{ T, W, input[0..packed_len(T, W)], reference, output });
                    return packed_len(T, W);
                }
            }
            unreachable;
        }

        pub fn dyn_undelta_pack(
            comptime T: type,
            noalias input: []align(1) const T,
            noalias base: *align(1) const [n_lanes(T)]T,
            noalias output: *[1024]T,
            width: usize,
        ) usize {
            inline for (0..@bitSizeOf(T) + 1) |W| {
                if (W == width) {
                    @call(.never_inline, undelta_pack, .{ T, W, input[0..packed_len(T, W)], base, output });
                    return packed_len(T, W);
                }
            }
            unreachable;
        }

fn packed_len(comptime T: type, comptime W: comptime_int) comptime_int {
    return 1024 * W / @bitSizeOf(T);
}

                pub fn bit_pack(
                    comptime T: type,
                    comptime W: comptime_int,
                    noalias input: *const [1024]T,
                    noalias output: *align(1) [packed_len(T, W)]T,
                ) void {
                    const N_BITS = @bitSizeOf(T);
                    const N_LANES = n_lanes(T);

                    for (0..n_lanes(T)) |lane| {
                        if (W == 0) {
                            return;
                        } else if (W == @bitSizeOf(T)) {
                            inline for (0..@bitSizeOf(T)) |row| {
                                const idx = index(row, lane);
                                output[N_LANES * row + lane] = input[idx];
                            }
                            return;
                        } else {
                            const mask_ = mask(T, W);

                            var tmp: T = 0;

                            inline for (0..N_BITS) |row| {
                                const idx = index(row, lane);
                                const src = input[idx] & mask_;

                                if (row == 0) {
                                    tmp = src;
                                } else {
                                    tmp |= src << (row * W) % N_BITS;
                                }

                                const curr_word: usize = (row * W) / N_BITS;
                                const next_word: usize = ((row + 1) * W) / N_BITS;

                                if (next_word > curr_word) {
                                    output[N_LANES * curr_word + lane] = tmp;
                                    const remaining_bits: usize = ((row + 1) * W) % N_BITS;
                                    tmp = src >> W - remaining_bits;
                                }
                            }
                        }
                    }
                }

                pub fn bit_unpack(
                    comptime T: type,
                    comptime W: comptime_int,
                    noalias input: *align(1) const [packed_len(T, W)]T,
                    noalias output: *[1024]T,
                ) void {
                    const N_BITS = @bitSizeOf(T);
                    const N_LANES = n_lanes(T);

                    for (0..N_LANES) |lane| {
                        if (W == 0) {
                            inline for (0..N_BITS) |row| {
                                const idx = index(row, lane);
                                output[idx] = 0;
                            }
                        } else if (W == N_BITS) {
                            inline for (0..N_BITS) |row| {
                                const idx = index(row, lane);
                                const src = input[N_LANES * row + lane];
                                output[idx] = src;
                            }
                        } else {
                            var src: T = input[lane];
                            var tmp: T = 0;

                            inline for (0..N_BITS) |row| {
                                const curr_word: usize = (row * W) / N_BITS;
                                const next_word: usize = ((row + 1) * W) / N_BITS;

                                const shift = (row * W) % N_BITS;

                                if (next_word > curr_word) {
                                    const remaining_bits = ((row + 1) * W) % N_BITS;
                                    const current_bits = W - remaining_bits;
                                    tmp = (src >> shift) & mask(T, current_bits);

                                    if (next_word < W) {
                                        src = input[N_LANES * next_word + lane];
                                        tmp |= (src & mask(T, remaining_bits)) << current_bits;
                                    }
                                } else {
                                    tmp = (src >> shift) & mask(T, W);
                                }

                                const idx = index(row, lane);
                                output[idx] = tmp;
                            }
                        }
                    }
                }

                pub fn for_pack(
                    comptime T: type,
                    comptime W: comptime_int,
                    noalias input: *const [1024]T,
                    reference: T,
                    noalias output: *align(1) [packed_len(T, W)]T,
                ) void {
                    const N_BITS = @bitSizeOf(T);
                    const N_LANES = n_lanes(T);

                    for (0..N_LANES) |lane| {
                        if (W == 0) {
                            return;
                        } else if (W == N_BITS) {
                            inline for (0..N_BITS) |row| {
                                const idx = index(row, lane);
                                output[N_LANES * row + lane] = input[idx] -% reference;
                            }
                            return;
                        } else {
                            const mask_ = mask(T, W);

                            var tmp: T = 0;

                            inline for (0..N_BITS) |row| {
                                const idx = index(row, lane);
                                const src = (input[idx] -% reference) & mask_;

                                if (row == 0) {
                                    tmp = src;
                                } else {
                                    tmp |= src << (row * W) % N_BITS;
                                }

                                const curr_word: usize = (row * W) / N_BITS;
                                const next_word: usize = ((row + 1) * W) / N_BITS;

                                if (next_word > curr_word) {
                                    output[N_LANES * curr_word + lane] = tmp;
                                    const remaining_bits: usize = ((row + 1) * W) % N_BITS;
                                    tmp = src >> W - remaining_bits;
                                }
                            }
                        }
                    }
                }

                pub fn for_unpack(
                    comptime T: type,
                    comptime W: comptime_int,
                    noalias input: *align(1) const [packed_len(T, W)]T,
                    reference: T,
                    noalias output: *[1024]T,
                ) void {
                    const N_BITS = @bitSizeOf(T);
                    const N_LANES = n_lanes(T);

                    for (0..N_LANES) |lane| {
                        if (W == 0) {
                            inline for (0..N_BITS) |row| {
                                const idx = index(row, lane);
                                output[idx] = reference;
                            }
                        } else if (W == N_BITS) {
                            inline for (0..N_BITS) |row| {
                                const idx = index(row, lane);
                                const src = input[N_LANES * row + lane];
                                output[idx] = src +% reference;
                            }
                        } else {
                            var src: T = input[lane];
                            var tmp: T = 0;

                            inline for (0..N_BITS) |row| {
                                const curr_word: usize = (row * W) / N_BITS;
                                const next_word: usize = ((row + 1) * W) / N_BITS;

                                const shift = (row * W) % N_BITS;

                                if (next_word > curr_word) {
                                    const remaining_bits = ((row + 1) * W) % N_BITS;
                                    const current_bits = W - remaining_bits;
                                    tmp = (src >> shift) & mask(T, current_bits);

                                    if (next_word < W) {
                                        src = input[N_LANES * next_word + lane];
                                        tmp |= (src & mask(T, remaining_bits)) << current_bits;
                                    }
                                } else {
                                    tmp = (src >> shift) & mask(T, W);
                                }

                                const idx = index(row, lane);
                                output[idx] = tmp +% reference;
                            }
                        }
                    }
                }

                pub fn undelta_pack(
                    comptime T: type,
                    comptime W: comptime_int,
                    noalias input: *align(1) const [packed_len(T, W)]T,
                    noalias base: *align(1) const [n_lanes(T)]T,
                    noalias output: *[1024]T,
                ) void {
                    const N_BITS = @bitSizeOf(T);
                    const N_LANES = n_lanes(T);
                    for (0..N_LANES) |lane| {
                        var prev: T = base[lane];

                        if (W == 0) {
                            inline for (0..N_BITS) |row| {
                                const idx = index(row, lane);
                                output[idx] = prev;
                            }
                        } else if (W == N_BITS) {
                            inline for (0..N_BITS) |row| {
                                const idx = index(row, lane);
                                const src = input[N_LANES * row + lane];

                                const next = src +% prev;
                                output[idx] = next;
                                prev = next;
                            }
                        } else {
                            var src: T = input[lane];
                            var tmp: T = 0;

                            inline for (0..N_BITS) |row| {
                                const curr_word: usize = (row * W) / N_BITS;
                                const next_word: usize = ((row + 1) * W) / N_BITS;

                                const shift = (row * W) % N_BITS;

                                if (next_word > curr_word) {
                                    const remaining_bits = ((row + 1) * W) % N_BITS;
                                    const current_bits = W - remaining_bits;
                                    tmp = (src >> shift) & mask(T, current_bits);

                                    if (next_word < W) {
                                        src = input[N_LANES * next_word + lane];
                                        tmp |= (src & mask(T, remaining_bits)) << current_bits;
                                    }
                                } else {
                                    tmp = (src >> shift) & mask(T, W);
                                }

                                const idx = index(row, lane);

                                const next = tmp +% prev;
                                output[idx] = next;
                                prev = next;
                            }
                        }
                    }
                }

// fn Test(comptime T: type) type {
//     return struct {
//         const std = @import("std");

//         const FL = FastLanes(T);

//         fn needed_width(noalias data: *const [1024]T) u7 {
//             var m = data[0];
//             for (0..1024) |idx| {
//                 m = @max(data[idx], m);
//             }
//             return @sizeOf(T) * 8 - @clz(m);
//         }

//         fn read_input(input: []const u8) ?[1024]T {
//             if (input.len < 1 + @sizeOf(T) * 1024) return null;

//             const has_zeroes = input[0] % 2 == 0;
//             var in: [1024]T = @bitCast(input[1 .. 1 + @sizeOf(T) * 1024].*);

//             if (!has_zeroes) {
//                 const replacement = in[0] +| 1;
//                 for (0..1024) |i| {
//                     if (in[i] == 0) {
//                         in[i] = replacement;
//                     }
//                 }
//             }

//             return in;
//         }

//         fn fuzz_rle(_: void, input: []const u8) anyerror!void {
//             const in = read_input(input) orelse return;

//             var rle_vals = std.mem.zeroes([1024]T);
//             var rle_idxs = std.mem.zeroes([1024]u16);

//             const len = FL.rle_encode(&in, &rle_vals, &rle_idxs);

//             var out = std.mem.zeroes([1024]T);

//             FL.rle_decode(rle_vals[0..len], &rle_idxs, &out);

//             try std.testing.expect(std.mem.eql(T, &out, &in));
//         }

//         test "fuzz_rle" {
//             try std.testing.fuzz({}, fuzz_rle, .{});
//         }

//         fn fuzz_bit_pack(_: void, input: []const u8) anyerror!void {
//             const in = read_input(input) orelse return;
//             const width = needed_width(&in);
//             var p = std.mem.zeroes([1024]T);
//             const packed_len = FL.dyn_bit_pack(&in, &p, width);

//             var out = std.mem.zeroes([1024]T);

//             const pl = FL.dyn_bit_unpack(&p, &out, width);

//             try std.testing.expectEqual(pl, packed_len);

//             try std.testing.expect(std.mem.eql(T, &out, &in));
//         }

//         test "fuzz_bit_pack" {
//             try std.testing.fuzz({}, fuzz_bit_pack, .{});
//         }

//         fn fuzz_delta_pack(_: void, input: []const u8) anyerror!void {
//             const in = read_input(input) orelse return;

//             var transposed = std.mem.zeroes([1024]T);
//             FL.transpose(&in, &transposed);

//             var untransposed = std.mem.zeroes([1024]T);
//             FL.untranspose(&transposed, &untransposed);

//             try std.testing.expect(std.mem.eql(T, &in, &untransposed));

//             const bases = transposed[0..FL.N_LANES];

//             var delta = std.mem.zeroes([1024]T);
//             FL.delta(&transposed, bases, &delta);

//             var undelta = std.mem.zeroes([1024]T);
//             FL.undelta(&delta, bases, &undelta);

//             try std.testing.expect(std.mem.eql(T, &undelta, &transposed));

//             const width = needed_width(&delta);

//             var p = std.mem.zeroes([1024]T);
//             const packed_len = FL.dyn_bit_pack(&delta, &p, width);

//             var unpacked = std.mem.zeroes([1024]T);
//             const pl = FL.dyn_bit_unpack(&p, &unpacked, width);

//             try std.testing.expectEqual(pl, packed_len);

//             try std.testing.expect(std.mem.eql(T, &delta, &unpacked));

//             var undelta_packed = std.mem.zeroes([1024]T);
//             const dpl = FL.dyn_undelta_pack(&p, bases, &undelta_packed, width);
//             try std.testing.expectEqual(dpl, packed_len);

//             try std.testing.expect(std.mem.eql(T, &undelta_packed, &transposed));
//         }

//         test "fuzz_delta_pack" {
//             try std.testing.fuzz({}, fuzz_delta_pack, .{});
//         }

//         fn fuzz_for_pack(_: void, input: []const u8) anyerror!void {
//             const in = read_input(input) orelse return;

//             const min, const max = std.mem.minMax(T, &in);
//             const width = (@sizeOf(T) * 8) - @clz(max - min);

//             var p = std.mem.zeroes([1024]T);
//             const packed_len = FL.dyn_for_pack(&in, min, &p, width);

//             var out = std.mem.zeroes([1024]T);

//             const pl = FL.dyn_for_unpack(&p, min, &out, width);

//             try std.testing.expectEqual(pl, packed_len);

//             try std.testing.expect(std.mem.eql(T, &out, &in));
//         }

//         test "fuzz_for_pack" {
//             try std.testing.fuzz({}, fuzz_for_pack, .{});
//         }
//     };
// }

// test "u8" {
//     _ = Test(u8);
// }

// test "u16" {
//     _ = Test(u16);
// }

// test "u32" {
//     _ = Test(u32);
// }

// test "u64" {
//     _ = Test(u64);
// }
