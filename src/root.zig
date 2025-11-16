const std = @import("std");

pub const Error = error {
    InvalidInput,
};

const BLOCK_LEN = 32;

fn pack_block(comptime T: type, block: [BLOCK_LEN]T, out: []u8) usize {
    const t_info = @typeInfo(T).int;

    if (t_info.bits > 256) {
        @compileError("only integers up to 256 bits are supported");
    }
 
    std.debug.assert(out.len >= @sizeOf(u8) + @sizeOf([BLOCK_LEN]T));
 
    const t_list = if (@typeInfo(T).int.signedness == .unsigned)
        .{u8, u16, u32, u64, u128, u256}
    else
        .{i8, i16, i32, i64, i128, i256};

    var max: T = block[0];
    for (0..BLOCK_LEN) |idx| {
        max = @max(max, block[idx]);
    }

    inline for (t_list) |out_t| {
        if (max <= std.math.maxInt(out_t)) {
            out[0] = @sizeOf(out_t);
            const out_p: [][@sizeOf(out_t)]u8 = @ptrCast(out[1..]);
            for (0..BLOCK_LEN) |idx| {
                out_p[idx] = @bitCast(@as(out_t, @truncate(block[idx])));
            }
            return @sizeOf([BLOCK_LEN]out_t);
        }
    }

    unreachable;
}

fn delta_pack_block(comptime T: type, block: [BLOCK_LEN]T, out: []u8) usize {
    std.debug.assert(out.len >= @sizeOf(T) + @sizeOf(u8) + @sizeOf([BLOCK_LEN]T));

    var delta: [BLOCK_LEN]T = undefined;

    const base = block[0];
    for (0..BLOCK_LEN) |idx| {
        delta[idx] = block[idx] -% base;
    }

    @as([][@sizeOf(T)]u8, @ptrCast(out))[0] = @bitCast(base);

    return @sizeOf(T) + pack_block(T, delta, out[@sizeOf(T)..]);
}

/// Delta encode and byte pack a list of integers.
///
/// `output.len` should be at least `delta_byte_pack_bound(T, input.len)`
///
/// returns the number of bytes written to the output.
pub fn delta_byte_pack(comptime T: type, noalias input: []const T, noalias output: []u8) Error!usize {
    if (input.len == 0) {
        return 0;
    }

    if (output.len < delta_byte_pack_bound(T, input.len)) {
        return Error.InvalidInput;
    }

    const num_blocks = input.len / BLOCK_LEN;

    const blocks = @as([]const [BLOCK_LEN]T, @ptrCast(input[0..num_blocks * BLOCK_LEN]));

    var offset: usize = 0;
    for (blocks) |block| {
        offset += delta_pack_block(T, block, output[offset..]);
    }

    var last_block: [BLOCK_LEN]T = std.mem.zeroes([BLOCK_LEN]T);
    const blocks_len = BLOCK_LEN * num_blocks;
    for (input[blocks_len..], 0..) |v, idx| {
        last_block[idx] = v;
    }
    offset += delta_pack_block(T, last_block, output[offset..]);

    return offset;
}

fn read_int(comptime T: type, input: []const u8) T {
    return @bitCast(@as([]const [@sizeOf(T)]u8, @ptrCast(input))[0]);
}

fn write_int(comptime T: type, val: T, output: []u8) void {
    @as([][@sizeOf(T)]u8, @ptrCast(output))[0] = @bitCast(val);
}

fn unpack_block(comptime PackT: type, comptime T: type, base: T, data: []const u8) [BLOCK_LEN]T {
    const Bits = [@sizeOf(PackT)]u8;

    const d: [BLOCK_LEN]PackT = @bitCast(@as([]const[BLOCK_LEN]Bits, @ptrCast(data))[0]);

    var out: [BLOCK_LEN]T = undefined;
    for (0..BLOCK_LEN) |idx| {
        out[idx] = d[idx] +% base;
    }

    return out;
}

fn load_block(comptime T: type, noalias offset: *usize, noalias input: []const u8) Error![BLOCK_LEN]T {
    const in = input[offset.*..];

    if (in.len < @sizeOf(T) + @sizeOf(u8)) {
        return Error.InvalidInput;
    }

    const base = read_int(T, input);
    const byte_width = read_int(u8, input[@sizeOf(T)..]);

    const data = in[@sizeOf(T) + @sizeOf(u8)..];

    if (data.len < byte_width * BLOCK_LEN) {
        return Error.InvalidInput;
    }

    offset += @sizeOf(T) + @sizeOf(u8) + byte_width * BLOCK_LEN;

    inline for (.{1, 2, 4, 8, 16}) |bits| {
        const pack_t = std.meta.Int(@typeInfo(T).int.signedness, bits);
        if (bits * 8 == byte_width) {
            return unpack_block(pack_t, T, base, data);
        }
    }

    return Error.InvalidInput;

}

pub fn delta_byte_unpack(comptime T: type, noalias input: []const u8, noalias output: []T) void {
    const num_blocks = output.len / BLOCK_LEN;

    var offset: usize = 0;

    const out_blocks: [][BLOCK_LEN]T = @ptrCast(output[0..num_blocks * BLOCK_LEN]);

    for (out_blocks) |*out_block| {
        out_block.* = try load_block(T, &offset, input);
    }

    const n_extra = output.len % BLOCK_LEN;

    if (n_extra > 0) {
        const extra_block = try load_block(T, &offset, input);

        const out = output[num_blocks * BLOCK_LEN..];

        for (0..n_extra) |idx| {
            out[idx] = extra_block[idx];
        }
    }
}

/// Calculate the maximum output size of `delta_byte_pack`
pub fn delta_byte_pack_bound(comptime T: type, len: usize) usize {
    const block_size = @sizeOf(T) + @sizeOf(u8) + @sizeOf([BLOCK_LEN]T);
    const num_blocks = (len + BLOCK_LEN - 1) / BLOCK_LEN;
    return block_size * num_blocks;
}

pub fn byte_pack() void {}

pub fn byte_unpack() void {}

pub fn byte_pack_bound() usize {}
