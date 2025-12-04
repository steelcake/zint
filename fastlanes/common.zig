const std = @import("std");

pub const FL_ORDER = [_]usize{ 0, 4, 2, 6, 1, 5, 3, 7 };

pub fn index(row: usize, lane: usize) usize {
    const o = row / 8;
    const s = row % 8;
    return (FL_ORDER[o] * 16) + (s * 128) + lane;
}

fn mask(comptime T: type, width: usize) T {
    return (1 << width) - 1;
}

pub fn packed_len(comptime T: type, comptime W: comptime_int) comptime_int {
    const n_bits = @sizeOf(T) * 8;
    return 1024 / W * n_bits;
}

pub inline fn pack(
    comptime T: type,
    comptime W: comptime_int,
    ctx: anytype,
    comptime kernel: fn (@TypeOf(ctx), idx: usize) T,
    noalias output: *[packed_len(T, W)]T,
    lane: usize,
) void {
    const n_bits = @sizeOf(T) * 8;
    const n_lanes = 1024 / n_bits;

    if (W == 0) {
        return;
    } else if (W == n_bits) {
        inline for (0..n_bits) |row| {
            const idx = index(row, lane);
            output[n_lanes * row + lane] = kernel(ctx, idx);
        }
        return;
    } else {
        const mask_ = mask(T, W);

        var tmp: T = 0;

        inline for (0..n_bits) |row| {
            const idx = index(row, lane);
            const src = kernel(ctx, idx) & mask_;

            if (row == 0) {
                tmp = src;
            } else {
                tmp |= src << (row * W) % n_bits;
            }

            const curr_word: usize = (row * W) / n_bits;
            const next_word: usize = ((row + 1) * W) / n_bits;

            if (next_word > curr_word) {
                output[n_lanes * curr_word + lane] = tmp;
                const remaining_bits: usize = ((row + 1) * W) % n_bits;
                tmp = src >> W - remaining_bits;
            }
        }
    }
}

pub inline fn unpack(
    comptime T: type,
    comptime W: comptime_int,
    ctx: anytype,
    comptime kernel: fn (@TypeOf(ctx), idx: usize, elem: T) void,
    noalias input: *const [packed_len(T, W)]T,
    lane: usize,
) void {
    const n_bits = @sizeOf(T) * 8;
    const n_lanes = 1024 / n_bits;

    if (W == 0) {
        inline for (0..n_bits) |row| {
            const idx = index(row, lane);
            kernel(ctx, idx, 0);
        }
    } else if (W == n_bits) {
        inline for (0..n_bits) |row| {
            const idx = index(row, lane);
            const src = input[n_lanes * row + lane];
            kernel(ctx, idx, src);
        }
    } else {
        var src: T = input[lane];
        var tmp: T = 0;

        inline for (0..n_bits) |row| {
            const curr_word: usize = (row * W) / n_bits;
            const next_word: usize = ((row + 1) * W) / n_bits;

            const shift = (row * W) % n_bits;

            if (next_word > curr_word) {
                const remaining_bits = ((row + 1) * W) % n_bits;
                const current_bits = W - remaining_bits;
                tmp = (src >> shift) & mask(current_bits);

                if (next_word < W) {
                    src = input[n_lanes * next_word + lane];
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
