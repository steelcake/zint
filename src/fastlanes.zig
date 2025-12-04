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
        const N_LANES = 1024 / N_BITS;

        fn mask(width: usize) T {
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
            noalias output: []u8,
            width: usize,
        ) usize {
            const out: []align(1) T = @ptrCast(output);
            inline for (0..N_BITS + 1) |W| {
                if (W == width) {
                    const P = Packer(W);
                    P.bit_pack(input, out[0..P.PACKED_LEN]);
                    return P.PACKED_LEN * @sizeOf(T);
                }
            }
            unreachable;
        }

        pub fn dyn_bit_unpack(
            noalias input: []const u8,
            noalias output: *[1024]T,
            width: usize,
        ) usize {
            const in: []align(1) const T = @ptrCast(input);
            inline for (0..N_BITS + 1) |W| {
                if (W == width) {
                    const P = Packer(W);
                    P.bit_unpack(in[0..P.PACKED_LEN], output);
                    return P.PACKED_LEN * @sizeOf(T);
                }
            }
            unreachable;
        }

        pub fn dyn_for_pack(
            noalias input: *const [1024]T,
            reference: T,
            noalias output: []u8,
            width: usize,
        ) usize {
            const out: []align(1) T = @ptrCast(output);
            inline for (0..N_BITS + 1) |W| {
                if (W == width) {
                    const P = Packer(W);
                    P.for_pack(input, reference, out[0..P.PACKED_LEN]);
                    return P.PACKED_LEN * @sizeOf(T);
                }
            }
            unreachable;
        }

        pub fn dyn_for_unpack(
            noalias input: []const u8,
            reference: T,
            noalias output: *[1024]T,
            width: usize,
        ) usize {
            const in: []align(1) const T = @ptrCast(input);
            inline for (0..N_BITS + 1) |W| {
                if (W == width) {
                    const P = Packer(W);
                    P.for_unpack(in[0..P.PACKED_LEN], reference, output);
                    return P.PACKED_LEN * @sizeOf(T);
                }
            }
            unreachable;
        }

        pub fn dyn_undelta_pack(
            noalias input: []const u8,
            noalias base: *const [N_LANES]T,
            noalias output: *[1024]T,
            width: usize,
        ) void {
            const in: []align(1) const T = @ptrCast(input);
            inline for (0..N_BITS + 1) |W| {
                if (W == width) {
                    const P = Packer(W);
                    P.undelta_pack(in[0..P.PACKED_LEN], base, output);
                    return P.PACKED_LEN * @sizeOf(T);
                }
            }
            unreachable;
        }

        pub fn Packer(comptime W: comptime_int) type {
            if (W > N_BITS) {
                @compileError("W can't be bigger than N_BITS");
            }

            return struct {
                const PACKED_LEN = 1024 / W * N_BITS;

                inline fn pack(
                    ctx: anytype,
                    comptime kernel: fn (@TypeOf(ctx), idx: usize) T,
                    noalias output: *[PACKED_LEN]T,
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
                        const mask_ = mask(T, W);

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
                    noalias input: *const [PACKED_LEN]T,
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
                    noalias output: *[PACKED_LEN]T,
                ) void {
                    const Kernel = struct {
                        fn kernel(noalias in: *const [1024]T, idx: usize) T {
                            return in[idx];
                        }
                    };

                    for (0..N_LANES) |lane| {
                        pack(T, W, input, Kernel.kernel, output, lane);
                    }
                }

                pub fn bit_unpack(
                    noalias input: *const [PACKED_LEN]T,
                    noalias output: *[1024]T,
                ) void {
                    const Kernel = struct {
                        fn kernel(noalias out: *[1024]T, idx: usize, elem: T) void {
                            out[idx] = elem;
                        }
                    };

                    for (0..N_LANES) |lane| {
                        unpack(T, W, output, Kernel.kernel, input, lane);
                    }
                }

                pub fn for_pack(
                    noalias input: *const [1024]T,
                    reference: T,
                    noalias output: *[PACKED_LEN]T,
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
                        pack(T, W, ctx, Kernel.kernel, output, lane);
                    }
                }

                pub fn for_unpack(
                    noalias input: *const [PACKED_LEN]T,
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
                        unpack(T, W, ctx, Kernel.kernel, input, lane);
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

                        unpack(T, W, ctx, Kernel.kernel, input, lane);
                    }
                }
            };
        }
    };
}
