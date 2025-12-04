const common = @import("common.zig");

pub fn for_pack(
    comptime T: type,
    comptime W: comptime_int,
    noalias input: *const [1024]T,
    reference: T,
    noalias output: *[common.packed_len(T, W)]T,
) void {
    const n_bits = @sizeOf(T) * 8;
    const n_lanes = 1024 / n_bits;

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

    for (0..n_lanes) |lane| {
        common.pack(T, W, ctx, Kernel.kernel, output, lane);
    }
}

pub fn for_unpack(
    comptime T: type,
    comptime W: comptime_int,
    noalias input: *const [common.packed_len(T, W)]T,
    reference: T,
    noalias output: *[1024]T,
) void {
    const n_bits = @sizeOf(T) * 8;
    const n_lanes = 1024 / n_bits;

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

    for (0..n_lanes) |lane| {
        common.unpack(T, W, ctx, Kernel.kernel, input, lane);
    }
}
