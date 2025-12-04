const common = @import("common.zig");

pub fn pack(
    comptime T: type,
    comptime W: comptime_int,
    noalias input: *const [1024]T,
    noalias output: *[common.packed_len(T, W)]T,
) void {
    const n_bits = @sizeOf(T) * 8;
    const n_lanes = 1024 / n_bits;

    const Kernel = struct {
        fn kernel(noalias in: *const [1024]T, idx: usize) T {
            return in[idx];
        }
    };

    for (0..n_lanes) |lane| {
        common.pack(T, W, input, Kernel.kernel, output, lane);
    }
}

pub fn unpack(
    comptime T: type,
    comptime W: comptime_int,
    noalias input: *const [common.packed_len(T, W)]T,
    noalias output: *[1024]T,
) void {
    const n_bits = @sizeOf(T) * 8;
    const n_lanes = 1024 / n_bits;

    const Kernel = struct {
        fn kernel(noalias out: *[1024]T, idx: usize, elem: T) void {
            out[idx] = elem;
        }
    };

    for (0..n_lanes) |lane| {
        common.unpack(T, W, output, Kernel.kernel, input, lane);
    }
}
