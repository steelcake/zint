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
            const n_blocks = len / 1024;
            const n_remainder = len % 1024;
            const max_block_size = @sizeOf(T) * 1024;

            // Layout of the output is:
            // - input length
            // - bit_width(u8) per block
            // - remainder values bit_width
            // - packed full blocks
            // - packed remainder values

            return @sizeOf(u32) + n_blocks + n_remainder * @sizeOf(T) + n_blocks * max_block_size;
        }

        fn needed_width(noalias range: T) u8 {
            return @sizeOf(T) * 8 - @clz(range);
        }

        pub fn bitpack_compress(noalias input: []const T, noalias output: []u8) Error!void {
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
            const byte_offset = 4 + n_whole_blocks + @intFromBool(n_remainder > 0);
            const out_t_len = 1024 * n_whole_blocks + n_remainder;
            const out: []align(1) T = @ptrCast(output[byte_offset .. byte_offset + out_t_len * @sizeOf(T)]);

            var offset: usize = 0;

            for (0..n_whole_blocks) |block_idx| {
                const start = block_idx * 1024;
                const block = input[start .. start + 1024];
            }

            // write remainder data
            if (n_remainder > 0) {
                if (n_whole_blocks == 0) {}

                const width = needed_width(std.mem.max(T, input[0..n_remainder]));
            }

            var offset: usize = 0;

            const whole_blocks: []const [1024]T = @ptrCast(input[0 .. n_whole_blocks * 1024]);
            for (whole_blocks) |block| {
                offset += FL.dyn_bit_pack(block, out[offset..]);
            }
        }

        pub fn bitpack_decompress() void {}

        pub fn forpack_compress() void {}
        pub fn forpack_decompress() void {}
        pub fn forpack_compress_bound() void {}

        pub fn delta_compress() void {}
        pub fn delta_decompress() void {}
        pub fn delta_compress_bound() void {}
    };
}
