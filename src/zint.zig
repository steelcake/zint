const std = @import("std");

const FastLanes = @import("./fastlanes.zig").FastLanes;
const ZigZag = @import("./zigzag.zig").ZigZag;
const ScalarBitpack = @import("./scalar_bitpack.zig").ScalarBitpack;

// pub fn Zint(comptime T: type) type {
//     _ = T;
//     return struct {};
// }

// fn Signed128(comptime T: type) type {}

// fn Unsigned128(comptime T: type) type {}

// fn Unsigned256(comptime T: type) type {}

// fn Signed(comptime T: type) type {
//     switch (T) {
//         i8, i16, i32, i64 => {},
//         else => @compileError("unexpected type"),
//     }

//     // const ZZ = ZigZag(T);

//     return struct {};
// }

pub const Error = error{InvalidInput};

fn Unsigned(comptime T: type) type {
    switch (T) {
        u8, u16, u32, u64 => {},
        else => @compileError("unexpected type."),
    }

    const FL = FastLanes(T);

    return struct {
        pub fn bitpack_compress_bound(len: u32) u32 {
            const n_blocks = (len + 1023) / 1024;
            const max_block_size = @sizeOf(T) * 1024;

            // Layout of the output is:
            // - input length
            // - remainder values bit_width(u8)
            // - bit_width(u8) per block
            // - packed remainder values
            // - packed full blocks
            //
            // We round the block count up to make sure we can comfortably write remainder data
            // to the output buffer without worrying of going out of bounds when the length
            // is < 1024.
            //
            // So we treat remainder data size as a full block of data size.

            return @sizeOf(u32) + n_blocks + max_block_size * n_blocks;
        }

        fn needed_width(noalias range: T) u8 {
            return @sizeOf(T) * 8 - @clz(range);
        }

        pub fn bitpack_compress(noalias input: []const T, noalias output: []u8) Error!usize {
            if (input.len > std.math.intMax(u32)) {
                return Error.InvalidInput;
            }
            const len: u32 = @intCast(input.len);

            const n_whole_blocks = len / 1024;
            const n_remainder = len % 1024;

            const output_bound = bitpack_compress_bound(len);
            if (output.len != output_bound) {
                return Error.InvalidInput;
            }

            // write input.length
            @as(*align(1) u32, @ptrCast(output[0..4])).* = len;

            // Start writing compressed data, skip some bytes for saving bit_widths later on
            // 4 is for the length, 1 is for the byte_width of remainder data, rest is for
            // byte widths of whole blocks.
            const byte_offset = 4 + 1 + n_whole_blocks;

            const out_t_len = 1024 * n_whole_blocks + n_remainder;
            const out: []align(1) T = @ptrCast(output[byte_offset .. byte_offset + out_t_len * @sizeOf(T)]);

            // number of T written so far, NOT number of bytes
            var offset: usize = 0;

            // write remainder data
            if (n_remainder > 0) {
                const max = std.mem.max(T, input[0..n_remainder]);
                const width = needed_width(max);

                output[4] = width;

                const remainder_packed_len = (@as(u64, width) * n_remainder + @sizeOf(T) - 1) / @sizeOf(T);

                if (n_whole_blocks > 0) {
                    // There is data after the remainder so we can treat it as a full block
                    // but then cut-off some at the end so the rest of the compression writes over the
                    // unused part
                    const block: *const [1024]T = input[0..1024];

                    _ = FL.dyn_bit_pack(block, out[offset..], width);

                    offset += remainder_packed_len;
                } else {
                    // there is no data after the remainder so we should copy it somewhere else
                    // and merge with other values in order to make it a whole block

                    // TODO: maybe shouldn't use stack here, can use extra output allocation or a context parameter
                    // for extra buffer space
                    var block = std.meta.zeroes([1024]T);

                    @memcpy(block[0..n_remainder], input[0..n_remainder]);

                    _ = FL.dyn_bit_pack(block, out[offset..], width);

                    return byte_offset + @sizeOf(T) * remainder_packed_len;
                }
            } else {
                output[4] = 0;
            }

            // Write whole blocks
            for (0..n_whole_blocks) |block_idx| {
                const start = block_idx * 1024;
                const block: *const [1024]T = input[start .. start + 1024];

                const max = max1024(block);
                const width = needed_width(max);

                output[5 + block_idx] = width;

                offset += FL.dyn_bit_pack(block, out[offset..], width);
            }

            return byte_offset + @sizeOf(T) * offset;
        }

        pub fn bitpack_decompress(noalias input: []const u8, noalias output: []T) void {
            if (output.len > std.math.intMax(u32)) {
                return Error.InvalidInput;
            }
            const len: u32 = @intCast(output.len);

            const n_whole_blocks = len / 1024;
            const n_remainder = len % 1024;

            if (len < 4 + 1 + n_whole_blocks) {
                return Error.InvalidInput;
            }

            const out_len: u32 = @as(*align(1) u32, @ptrCast(input[0..4])).*;
            if (len != out_len) {
                return Error.InvalidInput;
            }

            const remainder_width = input[4];

            const block_widths = input[5 .. 5 + n_whole_blocks];

            const remainder_packed_len = (@as(u64, n_remainder) * @as(u64, remainder_width) + @sizeOf(T) - 1) / @sizeOf(T);
            var total_packed_len = remainder_packed_len;
            for (block_widths) |w| {
                total_packed_len += @as(u64, w) * 1024 / @sizeOf(T);
            }

            const data_section_byte_offset = 4 + 1 + n_whole_blocks;

            const data_section_bytes = input[data_section_byte_offset..];
            if (data_section_bytes.len != total_packed_len * @sizeOf(T)) {
                return Error.InvalidInput;
            }

            const data_section: []align(1) const T = @ptrCast(data_section_bytes);

            // number of T read so far, NOT number of bytes
            var offset = 0;

            // read remainder data
            if (n_remainder > 0) {
                // How big a block of packed data with remainder_width would be
                const block_needed_len = @as(u64, remainder_width) * 1024 / @sizeOf(T);

                if (data_section.len >= block_needed_len) {
                    // there is enough data after the remainder section
                    // we can just unpack a full block with remainder_width
                    // and then discard the unneeded part

                    const block: []align(1) const T = data_section[0..block_needed_len];

                    _ = FL.dyn_bit_unpack(block, output[0..1024], remainder_width);

                    offset += remainder_packed_len;
                } else {
                    // there is not enough data after the remainder section
                    // we have to copy packed remainder data somewhere to create a 1024
                    // length input for unpacking
                    //
                    // WARNING: there can be actual full blocks after remainder data
                    // even if thre is not enough data to have a pseudo full block
                    // with remainder data width. This is because the packed full blocks that come
                    // after the packed remainder data might have lower bit width than the remainder data.

                    var block = std.mem.zeroes([1024]T);

                    // copy from remainder packed data
                    @memcpy(block[0..remainder_packed_len], data_section);

                    const block_slice = block[0..block_needed_len];

                    var out = std.mem.zeroes([1024]T);
                    _ = FL.dyn_bit_unpack(block_slice, &out, remainder_width);
                    @memcpy(output[0..n_remainder], out[0..n_remainder]);

                    offset += remainder_packed_len;
                }
            }

            // Read whole blocks
            for (0..n_whole_blocks) |block_idx| {
                const width = input[5 + block_idx];

                const out_offset = block_idx * 1024 + n_remainder;

                offset += FL.dyn_bit_unpack(data_section[offset..], output[out_offset .. out_offset + 1024], width);
            }

            std.debug.assert(offset == total_packed_len);
            std.debug.assert(offset == data_section.len);

            return data_section_byte_offset + @sizeOf(T) * offset;
        }

        // pub fn forpack_compress() void {}
        // pub fn forpack_decompress() void {}
        // pub fn forpack_compress_bound() void {}

        // pub fn delta_compress() void {}
        // pub fn delta_decompress() void {}
        // pub fn delta_compress_bound() void {}

        fn max1024(input: *const [1024]T) T {
            var m = input[0];
            for (0..1024) |i| {
                m = @max(input[i], m);
            }
            return m;
        }
    };
}
