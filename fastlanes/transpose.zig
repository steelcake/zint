const common = @import("./common.zig");

pub fn transpose_idx(idx: usize) usize {
    const lane = idx % 16;
    const order = (idx / 16) % 8;
    const row = idx / 128;

    return (lane * 64) + (common.FL_ORDER[order] * 8) + row;
}

pub fn transpose(
    comptime T: type,
    noalias input: *const [1024]T,
    noalias output: *[1024]T,
) void {
    for (0..1024) |i| {
        output[i] = input[transpose_idx(i)];
    }
}

pub fn untranspose(
    comptime T: type,
    noalias input: *const [1024]T,
    noalias output: *[1024]T,
) void {
    for (0..1024) |i| {
        output[transpose_idx(i)] = input[i];
    }
}
