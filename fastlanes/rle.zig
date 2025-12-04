pub fn encode(
    comptime T: type,
    noalias input: *const [1024]T,
    noalias rle_vals: *[1024]T,
    noalias rle_idxs: *[1024]u16,
) usize {
    var prev_val = input[0];
    rle_vals[0] = prev_val;
    rle_idxs[0] = 0;

    var pos_val: u16 = 0;
    var rle_val_idx: usize = 1;

    for (1..1024) |i| {
        const cur_val = input[i];

        if (cur_val != prev_val) {
            rle_vals[rle_val_idx] = cur_val;
            rle_val_idx += 1;
            pos_val += 1;
            prev_val = cur_val;
        }

        rle_idxs[i] = pos_val;
    }

    return rle_val_idx;
}

pub fn decode(
    comptime T: type,
    comptime I: type,
    noalias rle_vals: *const []T,
    noalias rle_idxs: *const [1024]I,
    noalias output: *[1024]T,
) void {
    for (0..1024) |i| {
        const idx = rle_idxs[i];
        output[i] = rle_vals[idx];
    }
}
