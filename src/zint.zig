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
    const SPack = ScalarBitpack(T);

    const N_BYTES = @sizeOf(T);
    const N_BITS = N_BYTES * 8;

    return struct {
        /// Returns the number of BYTES needed on the output buffer for compressing
        /// `len` elements. The compression will likely end up writing less data
        /// but this is needed to make sure the compression will succeed.
        pub fn bitpack_compress_bound(len: u32) u32 {
            const n_blocks = (len + 1023) / 1024;
            const max_block_size = N_BYTES * 1024;

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

        /// Compress the input integers into the output buffer.
        /// `output` should be at least `bitpack_compress_bound(input.len)` BYTES.
        ///
        /// Returns the number of bytes written to the output.
        pub fn bitpack_compress(noalias input: []const T, noalias output: []u8) Error!usize {
            if (input.len > std.math.maxInt(u32)) {
                return Error.InvalidInput;
            }
            const len: u32 = @intCast(input.len);

            const n_whole_blocks = len / 1024;
            const n_remainder = len % 1024;

            const output_bound = bitpack_compress_bound(len);
            if (output.len < output_bound) {
                return Error.InvalidInput;
            }

            // write input.length
            @as(*align(1) u32, @ptrCast(output[0..4])).* = len;

            // Start writing compressed data, skip some bytes for saving bit_widths later on
            // 4 is for the length, 1 is for the byte_width of remainder data, rest is for
            // byte widths of whole blocks.
            const byte_offset = 4 + 1 + n_whole_blocks;

            // We should have this much capacity even though we will write less
            const out_t_len = 1024 * (n_whole_blocks + 1);

            const out: []align(1) T = @ptrCast(output[byte_offset .. byte_offset + out_t_len * N_BYTES]);
            std.debug.assert(out.len == out_t_len);

            // number of T written so far, NOT number of bytes
            var offset: usize = 0;

            // write remainder data
            if (n_remainder > 0) {
                const max = std.mem.max(T, input[0..n_remainder]);
                const width = needed_width(max);

                output[4] = width;

                const remainder_packed_len = (@as(u64, width) * n_remainder + N_BITS - 1) / N_BITS;
                const n_written = SPack.bitpack(
                    input[0..n_remainder],
                    0,
                    out,
                    width,
                ) catch {
                    @panic("failed scalar pack, this should never happen");
                };
                std.debug.assert(n_written == remainder_packed_len);
                offset += remainder_packed_len;
            } else {
                output[4] = 0;
            }

            // Write whole blocks
            const input_blocks: []const [1024]T = @ptrCast(input[n_remainder..]);
            for (0..n_whole_blocks) |block_idx| {
                const block: *const [1024]T = &input_blocks[block_idx];

                const max = max1024(block);
                const width = needed_width(max);

                output[5 + block_idx] = width;

                offset += FL.dyn_bit_pack(block, out[offset..], width);
            }

            return byte_offset + N_BYTES * offset;
        }

        pub fn bitpack_decompress(noalias input: []const u8, noalias output: []T) Error!void {
            if (output.len > std.math.maxInt(u32)) {
                return Error.InvalidInput;
            }
            const len: u32 = @intCast(output.len);

            const n_whole_blocks = len / 1024;
            const n_remainder = len % 1024;

            if (input.len < 4 + 1 + n_whole_blocks) {
                return Error.InvalidInput;
            }

            const out_len: u32 = @as(*align(1) const u32, @ptrCast(input[0..4])).*;
            if (len != out_len) {
                return Error.InvalidInput;
            }

            const remainder_width = input[4];
            if (remainder_width > N_BITS) {
                return Error.InvalidInput;
            }

            const block_widths = input[5 .. 5 + n_whole_blocks];

            const remainder_packed_len = (@as(u64, n_remainder) * @as(u64, remainder_width) + N_BITS - 1) / N_BITS;
            var total_packed_len = remainder_packed_len;
            for (block_widths) |w| {
                if (w > N_BITS) {
                    return Error.InvalidInput;
                }
            }
            for (block_widths) |w| {
                total_packed_len += @as(u64, w) * 1024 / N_BITS;
            }

            const data_section_byte_offset = 4 + 1 + n_whole_blocks;

            const data_section_bytes = input[data_section_byte_offset..];
            if (data_section_bytes.len != total_packed_len * N_BYTES) {
                return Error.InvalidInput;
            }

            const data_section: []align(1) const T = @ptrCast(data_section_bytes);

            // number of T read so far, NOT number of bytes
            var offset: usize = 0;

            // read remainder data
            if (n_remainder > 0) {
                const n_read = try SPack.bitunpack(
                    data_section[0..remainder_packed_len],
                    0,
                    output[0..n_remainder],
                    remainder_width,
                );
                std.debug.assert(n_read == remainder_packed_len);
                offset += remainder_packed_len;
            }

            // Read whole blocks
            for (0..n_whole_blocks) |block_idx| {
                const width = input[5 + block_idx];

                const out_offset = output[block_idx * 1024 + n_remainder ..];

                offset += FL.dyn_bit_unpack(data_section[offset..], out_offset[0..1024], width);
            }

            std.debug.assert(offset == total_packed_len);
            std.debug.assert(offset == data_section.len);

            const end = data_section_byte_offset + @sizeOf(T) * offset;
            if (end != input.len) {
                return Error.InvalidInput;
            }

            return;
        }

        /// Returns the number of BYTES needed on the output buffer for compressing
        /// `len` elements. The compression will likely end up writing less data
        /// but this is needed to make sure the compression will succeed.
        pub fn forpack_compress_bound(len: u32) u32 {
            const n_blocks = (len + 1023) / 1024;
            const max_block_size = N_BYTES * 1024;

            // Layout of the output is:
            // - input length
            // - remainder values bit_width(u8)
            // - bit_width(u8) per block
            // - packed remainder values, a reference value before the packed values
            // - packed full blocks, a reference value before every block
            //
            // We round the block count up to make sure we can comfortably write remainder data
            // to the output buffer without worrying of going out of bounds when the length
            // is < 1024.
            //
            // So we treat remainder data size as a full block of data size.

            return @sizeOf(u32) + (N_BYTES + 1) * n_blocks + max_block_size * n_blocks;
        }

        /// Compress the input integers into the output buffer.
        /// `output` should be at least `forpack_compress_bound(input.len)` BYTES.
        ///
        /// Returns the number of bytes written to the output.
        pub fn forpack_compress(noalias input: []const T, noalias output: []u8) Error!usize {
            if (input.len > std.math.maxInt(u32)) {
                return Error.InvalidInput;
            }
            const len: u32 = @intCast(input.len);

            const n_whole_blocks = len / 1024;
            const n_remainder = len % 1024;

            const output_bound = forpack_compress_bound(len);
            if (output.len < output_bound) {
                return Error.InvalidInput;
            }

            // write input.length
            @as(*align(1) u32, @ptrCast(output[0..4])).* = len;

            // Start writing compressed data, skip some bytes for saving bit_widths later on
            // 4 is for the length, 1 is for the byte_width of remainder data, rest is for
            // byte widths of whole blocks.
            const byte_offset = 4 + 1 + n_whole_blocks;

            // We should have this much capacity even though we will write less
            const out_t_len = (1024 + 1) * (n_whole_blocks + 1);

            const out: []align(1) T = @ptrCast(output[byte_offset .. byte_offset + out_t_len * N_BYTES]);
            std.debug.assert(out.len == out_t_len);

            // number of T written so far, NOT number of bytes
            var offset: usize = 0;

            // write remainder data
            if (n_remainder > 0) {
                const min, const max = std.mem.minMax(T, input[0..n_remainder]);
                const width = needed_width(max - min);

                output[4] = width;

                out[0] = min;
                offset += 1;

                const remainder_packed_len = (@as(u64, width) * n_remainder + N_BITS - 1) / N_BITS;
                const n_written = SPack.bitpack(
                    input[0..n_remainder],
                    min,
                    out[offset..],
                    width,
                ) catch {
                    @panic("failed scalar pack, this should never happen");
                };
                std.debug.assert(n_written == remainder_packed_len);
                offset += remainder_packed_len;
            } else {
                out[offset] = 0;
                offset += 1;

                output[4] = 0;
            }

            // Write whole blocks
            const input_blocks: []const [1024]T = @ptrCast(input[n_remainder..]);
            for (0..n_whole_blocks) |block_idx| {
                const block: *const [1024]T = &input_blocks[block_idx];

                const min, const max = minmax1024(block);
                const width = needed_width(max - min);

                output[5 + block_idx] = width;

                out[offset] = min;
                offset += 1;

                offset += FL.dyn_for_pack(block, min, out[offset..], width);
            }

            return byte_offset + N_BYTES * offset;
        }

        pub fn forpack_decompress(noalias input: []const u8, noalias output: []T) Error!void {
            if (output.len > std.math.maxInt(u32)) {
                return Error.InvalidInput;
            }
            const len: u32 = @intCast(output.len);

            const n_whole_blocks = len / 1024;
            const n_remainder = len % 1024;

            if (input.len < 4 + 1 + n_whole_blocks) {
                return Error.InvalidInput;
            }

            const out_len: u32 = @as(*align(1) const u32, @ptrCast(input[0..4])).*;
            if (len != out_len) {
                return Error.InvalidInput;
            }

            const remainder_width = input[4];
            if (remainder_width > N_BITS) {
                return Error.InvalidInput;
            }

            const block_widths = input[5 .. 5 + n_whole_blocks];

            const remainder_packed_len = (@as(u64, n_remainder) * @as(u64, remainder_width) + N_BITS - 1) / N_BITS;
            var total_packed_len = remainder_packed_len + 1;
            for (block_widths) |w| {
                if (w > N_BITS) {
                    return Error.InvalidInput;
                }
            }
            for (block_widths) |w| {
                total_packed_len += @as(u64, w) * 1024 / N_BITS + 1;
            }

            const data_section_byte_offset = 4 + 1 + n_whole_blocks;

            const data_section_bytes = input[data_section_byte_offset..];
            if (data_section_bytes.len != total_packed_len * N_BYTES) {
                return Error.InvalidInput;
            }

            const data_section: []align(1) const T = @ptrCast(data_section_bytes);

            // number of T read so far, NOT number of bytes
            var offset: usize = 0;

            // read remainder data
            if (n_remainder > 0) {
                const ref = data_section[offset];
                offset += 1;

                const n_read = try SPack.bitunpack(
                    data_section[offset .. offset + remainder_packed_len],
                    ref,
                    output[0..n_remainder],
                    remainder_width,
                );
                std.debug.assert(n_read == remainder_packed_len);
                offset += remainder_packed_len;
            } else {
                // as if we read remainder data min value
                offset += 1;
            }

            // Read whole blocks
            for (0..n_whole_blocks) |block_idx| {
                const width = input[5 + block_idx];

                const out_offset = output[block_idx * 1024 + n_remainder ..];

                const ref = data_section[offset];
                offset += 1;

                offset += FL.dyn_for_unpack(data_section[offset..], ref, out_offset[0..1024], width);
            }

            std.debug.assert(offset == total_packed_len);
            std.debug.assert(offset == data_section.len);

            const end = data_section_byte_offset + @sizeOf(T) * offset;
            if (end != input.len) {
                return Error.InvalidInput;
            }

            return;
        }

        // pub fn delta_compress_bound() void {}
        // pub fn delta_compress() void {}
        // pub fn delta_decompress() void {}

        fn max1024(input: *const [1024]T) T {
            var m = input[0];
            for (0..1024) |i| {
                m = @max(input[i], m);
            }
            return m;
        }

        fn minmax1024(input: *const [1024]T) struct { T, T } {
            var min = input[0];
            var max = input[0];

            for (0..1024) |i| {
                min = @min(input[i], min);
                max = @max(input[i], max);
            }

            return .{ min, max };
        }

        fn needed_width(range: T) u8 {
            return N_BITS - @clz(range);
        }
    };
}

fn TestUnsigned(comptime T: type) type {
    return struct {
        const MAX_NUM_INTS = 123321;

        const U = Unsigned(T);

        const Context = struct {
            const page_allocator = std.heap.page_allocator;

            input: []T,
            compressed: []u8,
            output: []T,

            fn compress_bound() usize {
                return @max(
                    U.bitpack_compress_bound(MAX_NUM_INTS),
                    U.forpack_compress_bound(MAX_NUM_INTS),
                );
            }

            fn init() Context {
                const input = page_allocator.alloc(T, MAX_NUM_INTS) catch unreachable;
                const output = page_allocator.alloc(T, MAX_NUM_INTS) catch unreachable;
                const compressed = page_allocator.alloc(
                    u8,
                    compress_bound(),
                ) catch unreachable;

                @memset(input, 0);
                @memset(output, 0);
                @memset(compressed, 0);

                return .{
                    .input = input,
                    .output = output,
                    .compressed = compressed,
                };
            }

            fn deinit(self: *Context) void {
                page_allocator.free(self.input);
                page_allocator.free(self.output);
                page_allocator.free(self.compressed);

                self.input = &.{};
                self.output = &.{};
                self.compressed = &.{};
            }
        };

        fn read_input(typed_input: []T, input: []const u8) []T {
            std.debug.assert(typed_input.len == MAX_NUM_INTS);

            const input_num_ints = input.len / @sizeOf(T);
            const num_ints = @min(input_num_ints, MAX_NUM_INTS);
            const num_bytes = @sizeOf(T) * num_ints;
            @memcpy(@as([]u8, @ptrCast(typed_input))[0..num_bytes], input[0..num_bytes]);

            return typed_input[0..num_ints];
        }

        fn roundtrip_bitpack(ctx: Context, input: []const u8) anyerror!void {
            const in = read_input(ctx.input, input);

            const compressed_len = try U.bitpack_compress(
                in,
                ctx.compressed,
            );

            std.debug.assert(
                compressed_len <= U.bitpack_compress_bound(MAX_NUM_INTS),
            );

            U.bitpack_decompress(
                ctx.compressed[0..compressed_len],
                ctx.output[0..in.len],
            ) catch unreachable;

            try std.testing.expectEqualSlices(T, in, ctx.output[0..in.len]);
        }

        test roundtrip_bitpack {
            var ctx = Context.init();
            defer ctx.deinit();
            try std.testing.fuzz(ctx, roundtrip_bitpack, .{});
        }

        fn roundtrip_forpack(ctx: Context, input: []const u8) anyerror!void {
            const in = read_input(ctx.input, input);

            const compressed_len = try U.forpack_compress(
                in,
                ctx.compressed,
            );

            std.debug.assert(
                compressed_len <= U.forpack_compress_bound(MAX_NUM_INTS),
            );

            U.forpack_decompress(
                ctx.compressed[0..compressed_len],
                ctx.output[0..in.len],
            ) catch unreachable;

            try std.testing.expectEqualSlices(T, in, ctx.output[0..in.len]);
        }

        test roundtrip_forpack {
            var ctx = Context.init();
            defer ctx.deinit();
            try std.testing.fuzz(ctx, roundtrip_forpack, .{});
        }

        fn garbage_bitpack(ctx: []T, input: []const u8) anyerror!void {
            const output = ctx;

            if (input.len < 4) {
                U.bitpack_decompress(input, output) catch return;
            } else {
                const n_ints: u32 = @bitCast(@as([*]const [4]u8, @ptrCast(input.ptr))[0]);
                const num_ints = @min(@as(usize, n_ints), MAX_NUM_INTS);
                U.bitpack_decompress(input[4..], output[0..num_ints]) catch return;
            }
        }

        test garbage_bitpack {
            var ctx = Context.init();
            defer ctx.deinit();
            try std.testing.fuzz(ctx.input, garbage_bitpack, .{});
        }
    };
}

test Unsigned {
    _ = TestUnsigned(u8);
    _ = TestUnsigned(u16);
    _ = TestUnsigned(u32);
    _ = TestUnsigned(u64);
}
