const common = @import("./common.zig");

pub fn delta(
    comptime T: type,
    noalias input: *const [1024]T,
    noalias base: *const [1024 / (@sizeOf(T) * 8)]T,
    noalias output: *[1024]T,
) void {
    const n_bits = @sizeOf(T) * 8;
    const n_lanes = 1024 / n_bits;

    for (0..n_lanes) |lane| {
        var prev = base[lane];
        inline for (0..n_bits) |row| {
            const idx = common.index(row, lane);
            const next = input[idx];
            output[idx] = next -% prev;
            prev = next;
        }
    }
}

pub fn undelta(
    comptime T: type,
    noalias input: *const [1024]T,
    noalias base: *const [1024 / (@sizeOf(T) * 8)]T,
    noalias output: *[1024]T,
) void {
    const n_bits = @sizeOf(T) * 8;
    const n_lanes = 1024 / n_bits;

    for (0..n_lanes) |lane| {
        var prev = base[lane];
        inline for (0..n_bits) |row| {
            const idx = common.index(row, lane);
            const next = input[idx] +% prev;
            output[idx] = next;
            prev = next;
        }
    }
}

pub fn undelta_pack(
    comptime T: type,
    comptime W: comptime_int,
    noalias input: *const [common.packed_len(T, W)]T,
    noalias base: *const [1024 / (@sizeOf(T) * 8)]T,
    noalias output: *[1024]T,
) void {
    const n_bits = @sizeOf(T) * 8;
    const n_lanes = 1024 / n_bits;

    for (0..n_lanes) |lane| {
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

        common.unpack(T, W, ctx, Kernel.kernel, input, lane);
    }
}
