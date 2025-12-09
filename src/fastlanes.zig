//! This file implements bit packing and delta coding.
//!
//! Code is a port of fastlanes Rust implementation by spiraldb:
//! github.com/spiraldb/fastlanes
//!
//! It was specifically ported from commit hash 8b655cf.
//!
//! Same license as the original code can be found in project root `LICENSE-APACHE`

const FL_ORDER = [_]usize{ 0, 4, 2, 6, 1, 5, 3, 7 };

fn index(row: usize, lane: usize) usize {
    const o = row / 8;
    const s = row % 8;
    return (FL_ORDER[o] * 16) + (s * 128) + lane;
}

fn transpose_idx(idx: usize) usize {
    const lane = idx % 16;
    const order = (idx / 16) % 8;
    const row = idx / 128;

    return (lane * 64) + (FL_ORDER[order] * 8) + row;
}

pub fn FastLanes(comptime T: type) type {
    switch (T) {
        u8, u16, u32, u64 => {},
        else => @compileError("FastLanes only supports u8, u16, u32, u64."),
    }

    return struct {
        const N_BITS = @sizeOf(T) * 8;
        pub const N_LANES = 1024 / N_BITS;

        fn mask(width: comptime_int) T {
            return (1 << width) - 1;
        }

        pub fn delta(
            noalias input: *const [1024]T,
            noalias base: *const [N_LANES]T,
            noalias output: *[1024]T,
        ) void {
            for (0..N_LANES) |lane| {
                var prev = base[lane];
                inline for (0..N_BITS) |row| {
                    const idx = index(row, lane);
                    const next = input[idx];
                    output[idx] = next -% prev;
                    prev = next;
                }
            }
        }

        pub fn undelta(
            noalias input: *const [1024]T,
            noalias base: *const [N_LANES]T,
            noalias output: *[1024]T,
        ) void {
            for (0..N_LANES) |lane| {
                var prev = base[lane];
                inline for (0..N_BITS) |row| {
                    const idx = index(row, lane);
                    const next = input[idx] +% prev;
                    output[idx] = next;
                    prev = next;
                }
            }
        }

        pub fn transpose(
            noalias input: *const [1024]T,
            noalias output: *[1024]T,
        ) void {
            for (0..1024) |i| {
                output[i] = input[transpose_idx(i)];
            }
        }

        pub fn untranspose(
            noalias input: *const [1024]T,
            noalias output: *[1024]T,
        ) void {
            for (0..1024) |i| {
                output[transpose_idx(i)] = input[i];
            }
        }

        pub fn dyn_bit_pack(
            noalias input: *const [1024]T,
            noalias output: []align(1) T,
            width: usize,
        ) usize {
            inline for (0..N_BITS + 1) |W| {
                if (W == width) {
                    const P = Packer(W);
                    P.bit_pack(input, output[0..P.PACKED_LEN]);
                    return P.PACKED_LEN;
                }
            }
            unreachable;
        }

        pub fn dyn_bit_unpack(
            noalias input: []align(1) const T,
            noalias output: *[1024]T,
            width: usize,
        ) usize {
            inline for (0..N_BITS + 1) |W| {
                if (W == width) {
                    const P = Packer(W);
                    P.bit_unpack(input[0..P.PACKED_LEN], output);
                    return P.PACKED_LEN;
                }
            }
            unreachable;
        }

        pub fn dyn_for_pack(
            noalias input: *const [1024]T,
            reference: T,
            noalias output: []align(1) T,
            width: usize,
        ) usize {
            inline for (0..N_BITS + 1) |W| {
                if (W == width) {
                    const P = Packer(W);
                    P.for_pack(input, reference, output[0..P.PACKED_LEN]);
                    return P.PACKED_LEN;
                }
            }
            unreachable;
        }

        pub fn dyn_for_unpack(
            noalias input: []align(1) const T,
            reference: T,
            noalias output: *[1024]T,
            width: usize,
        ) usize {
            inline for (0..N_BITS + 1) |W| {
                if (W == width) {
                    const P = Packer(W);
                    P.for_unpack(input[0..P.PACKED_LEN], reference, output);
                    return P.PACKED_LEN;
                }
            }
            unreachable;
        }

        pub fn dyn_undelta_pack(
            noalias input: []const T,
            noalias base: *const [N_LANES]T,
            noalias output: *[1024]T,
            width: usize,
        ) usize {
            inline for (0..N_BITS + 1) |W| {
                if (W == width) {
                    const P = Packer(W);
                    P.undelta_pack(input[0..P.PACKED_LEN], base, output);
                    return P.PACKED_LEN;
                }
            }
            unreachable;
        }

        pub fn Packer(comptime W: comptime_int) type {
            if (W > N_BITS) {
                @compileError("W can't be bigger than N_BITS");
            }

            return struct {
                const PACKED_LEN = 1024 * W / N_BITS;

                inline fn pack(
                    ctx: anytype,
                    comptime kernel: fn (@TypeOf(ctx), idx: usize) T,
                    noalias output: *align(1) [PACKED_LEN]T,
                    lane: usize,
                ) void {
                    if (W == 0) {
                        return;
                    } else if (W == N_BITS) {
                        inline for (0..N_BITS) |row| {
                            const idx = index(row, lane);
                            output[N_LANES * row + lane] = kernel(ctx, idx);
                        }
                        return;
                    } else {
                        const mask_ = mask(W);

                        var tmp: T = 0;

                        inline for (0..N_BITS) |row| {
                            const idx = index(row, lane);
                            const src = kernel(ctx, idx) & mask_;

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

                inline fn unpack(
                    ctx: anytype,
                    comptime kernel: fn (@TypeOf(ctx), idx: usize, elem: T) void,
                    noalias input: *align(1) const [PACKED_LEN]T,
                    lane: usize,
                ) void {
                    if (W == 0) {
                        inline for (0..N_BITS) |row| {
                            const idx = index(row, lane);
                            kernel(ctx, idx, 0);
                        }
                    } else if (W == N_BITS) {
                        inline for (0..N_BITS) |row| {
                            const idx = index(row, lane);
                            const src = input[N_LANES * row + lane];
                            kernel(ctx, idx, src);
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
                                tmp = (src >> shift) & mask(current_bits);

                                if (next_word < W) {
                                    src = input[N_LANES * next_word + lane];
                                    tmp |= (src & mask(remaining_bits)) << current_bits;
                                }
                            } else {
                                tmp = (src >> shift) & mask(W);
                            }

                            const idx = index(row, lane);
                            kernel(ctx, idx, tmp);
                        }
                    }
                }

                pub fn bit_pack(
                    noalias input: *const [1024]T,
                    noalias output: *align(1) [PACKED_LEN]T,
                ) void {
                    const Kernel = struct {
                        fn kernel(in: *const [1024]T, idx: usize) T {
                            return in[idx];
                        }
                    };

                    for (0..N_LANES) |lane| {
                        pack(input, Kernel.kernel, output, lane);
                    }
                }

                pub fn bit_unpack(
                    noalias input: *align(1) const [PACKED_LEN]T,
                    noalias output: *[1024]T,
                ) void {
                    const Kernel = struct {
                        fn kernel(out: *[1024]T, idx: usize, elem: T) void {
                            out[idx] = elem;
                        }
                    };

                    for (0..N_LANES) |lane| {
                        unpack(output, Kernel.kernel, input, lane);
                    }
                }

                pub fn for_pack(
                    noalias input: *const [1024]T,
                    reference: T,
                    noalias output: *align(1) [PACKED_LEN]T,
                ) void {
                    const Ctx = struct {
                        ref: T,
                        in: *const [1024]T,
                    };

                    const ctx = Ctx{
                        .ref = reference,
                        .in = input,
                    };

                    const Kernel = struct {
                        fn kernel(c: Ctx, idx: usize) T {
                            return c.in[idx] -% c.ref;
                        }
                    };

                    for (0..N_LANES) |lane| {
                        pack(ctx, Kernel.kernel, output, lane);
                    }
                }

                pub fn for_unpack(
                    noalias input: *align(1) const [PACKED_LEN]T,
                    reference: T,
                    noalias output: *[1024]T,
                ) void {
                    const Ctx = struct {
                        ref: T,
                        out: *[1024]T,
                    };

                    const ctx = Ctx{
                        .ref = reference,
                        .out = output,
                    };

                    const Kernel = struct {
                        fn kernel(c: Ctx, idx: usize, elem: T) void {
                            c.out[idx] = elem +% c.ref;
                        }
                    };

                    for (0..N_LANES) |lane| {
                        unpack(ctx, Kernel.kernel, input, lane);
                    }
                }

                pub fn undelta_pack(
                    noalias input: *const [PACKED_LEN]T,
                    noalias base: *const [N_LANES]T,
                    noalias output: *[1024]T,
                ) void {
                    for (0..N_LANES) |lane| {
                        var prev: T = base[lane];

                        const Ctx = struct {
                            prev: *T,
                            out: *[1024]T,
                        };

                        const ctx = Ctx{
                            .prev = &prev,
                            .out = output,
                        };

                        const Kernel = struct {
                            fn kernel(c: Ctx, idx: usize, elem: T) void {
                                const next = elem +% c.prev.*;
                                c.out[idx] = next;
                                c.prev.* = next;
                            }
                        };

                        unpack(ctx, Kernel.kernel, input, lane);
                    }
                }
            };
        }
    };
}

fn Test(comptime T: type) type {
    return struct {
        const std = @import("std");

        const FL = FastLanes(T);

        fn needed_width(noalias data: *const [1024]T) u7 {
            var m = data[0];
            for (0..1024) |idx| {
                m = @max(data[idx], m);
            }
            return @sizeOf(T) * 8 - @clz(m);
        }

        fn read_input(input: []const u8) ?[1024]T {
            if (input.len < 1 + @sizeOf(T) * 1024) return null;

            const has_zeroes = input[0] % 2 == 0;
            var in: [1024]T = @bitCast(input[1 .. 1 + @sizeOf(T) * 1024].*);

            if (!has_zeroes) {
                const replacement = in[0] +| 1;
                for (0..1024) |i| {
                    if (in[i] == 0) {
                        in[i] = replacement;
                    }
                }
            }

            return in;
        }

        fn fuzz_bit_pack(_: void, input: []const u8) anyerror!void {
            const in = read_input(input) orelse return;
            const width = needed_width(&in);
            var p = std.mem.zeroes([1024]T);
            const packed_len = FL.dyn_bit_pack(&in, &p, width);

            var out = std.mem.zeroes([1024]T);

            const pl = FL.dyn_bit_unpack(&p, &out, width);

            try std.testing.expectEqual(pl, packed_len);

            try std.testing.expect(std.mem.eql(T, &out, &in));
        }

        test "fuzz_bit_pack" {
            try std.testing.fuzz({}, fuzz_bit_pack, .{});
        }

        fn fuzz_delta_pack(_: void, input: []const u8) anyerror!void {
            const in = read_input(input) orelse return;

            var transposed = std.mem.zeroes([1024]T);
            FL.transpose(&in, &transposed);

            var untransposed = std.mem.zeroes([1024]T);
            FL.untranspose(&transposed, &untransposed);

            try std.testing.expect(std.mem.eql(T, &in, &untransposed));

            const bases = transposed[0..FL.N_LANES];

            var delta = std.mem.zeroes([1024]T);
            FL.delta(&transposed, bases, &delta);

            var undelta = std.mem.zeroes([1024]T);
            FL.undelta(&delta, bases, &undelta);

            try std.testing.expect(std.mem.eql(T, &undelta, &transposed));

            const width = needed_width(&delta);

            var p = std.mem.zeroes([1024]T);
            const packed_len = FL.dyn_bit_pack(&delta, &p, width);

            var unpacked = std.mem.zeroes([1024]T);
            const pl = FL.dyn_bit_unpack(&p, &unpacked, width);

            try std.testing.expectEqual(pl, packed_len);

            try std.testing.expect(std.mem.eql(T, &delta, &unpacked));

            var undelta_packed = std.mem.zeroes([1024]T);
            const dpl = FL.dyn_undelta_pack(&p, bases, &undelta_packed, width);
            try std.testing.expectEqual(dpl, packed_len);

            try std.testing.expect(std.mem.eql(T, &undelta_packed, &transposed));
        }

        test "fuzz_delta_pack" {
            try std.testing.fuzz({}, fuzz_delta_pack, .{});
        }

        fn fuzz_for_pack(_: void, input: []const u8) anyerror!void {
            const in = read_input(input) orelse return;

            const min, const max = std.mem.minMax(T, &in);
            const width = (@sizeOf(T) * 8) - @clz(max - min);

            var p = std.mem.zeroes([1024]T);
            const packed_len = FL.dyn_for_pack(&in, min, &p, width);

            var out = std.mem.zeroes([1024]T);

            const pl = FL.dyn_for_unpack(&p, min, &out, width);

            try std.testing.expectEqual(pl, packed_len);

            try std.testing.expect(std.mem.eql(T, &out, &in));
        }

        test "fuzz_for_pack" {
            try std.testing.fuzz({}, fuzz_for_pack, .{});
        }
    };
}

test "u8" {
    _ = Test(u8);
}

test "u16" {
    _ = Test(u16);
}

test "u32" {
    _ = Test(u32);
}

test "u64" {
    _ = Test(u64);
}
