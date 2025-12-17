const std = @import("std");
const Allocator = std.mem.Allocator;

const FastLanes = @import("./fastlanes.zig").FastLanes;
const ZigZag = @import("./zigzag.zig").ZigZag;
const ScalarBitpack = @import("./scalar_bitpack.zig").ScalarBitpack;

const CONTEXT_SIZE = 1 << 16; // 64KB context is needed in case of delta compressing 256bit integers

pub const Error = error{InvalidInput};

/// Context that is used for compression and decompression operations.
///
/// This is a fairly large allocation (currently 64KB). So it is recommended to re-use this.
pub const Ctx = struct {
    buf: []align(64) u8,

    pub fn init(alloc: Allocator) error{OutOfMemory}!Ctx {
        const buf = try alloc.alignedAlloc(u8, std.mem.Alignment.fromByteUnits(64), CONTEXT_SIZE);

        return .{
            .buf = buf,
        };
    }

    pub fn deinit(self: Ctx, alloc: Allocator) void {
        alloc.free(self.buf);
    }

    fn typed(self: Ctx, comptime T: type) []align(64) T {
        return @ptrCast(self.buf);
    }
};

/// Structure to compress/decompress integers with type `T`
pub fn Zint(comptime T: type) type {
    return struct {
        /// Compression function will expect this amount of bytes of capacity on the
        /// output buffer.
        pub fn bitpack_compress_bound(len: u32) u32 {
            return switch (T) {
                i128, u128 => Impl128(T).bitpack_compress_bound(len),
                i256, u256 => Impl256(T).bitpack_compress_bound(len),
                else => Impl(T).bitpack_compress_bound(len),
            };
        }

        /// Compress integers from `input` to `output` buffer.
        ///
        /// Will error if `output.len < bitpack_compress_bound(input.len)` or if
        /// `input.len >= maxInt(u32)`.
        ///
        /// Returns the number of bytes written to the `output`.
        pub fn bitpack_compress(
            ctx: Ctx,
            noalias input: []const T,
            noalias output: []u8,
        ) Error!usize {
            const buf = ctx.typed(T);

            return switch (T) {
                i128, u128 => Impl128(T).bitpack_compress(
                    buf[0..1024],
                    buf[1024..1536],
                    input,
                    output,
                ),
                i256, u256 => Impl256(T).bitpack_compress(
                    buf[0..1024],
                    buf[1024..1280],
                    input,
                    output,
                ),
                else => Impl(T).bitpack_compress(
                    buf[0..1024],
                    input,
                    output,
                ),
            };
        }

        /// Decompress integers from `input` buffer to `output`.
        ///
        /// Expects enough data to decompress `output.len` integers from `input`.
        ///
        /// Returns the byte count of data read from `input` buffer.
        ///
        /// Will error if input data is invalid or short.
        pub fn bitpack_decompress(
            ctx: Ctx,
            noalias input: []const u8,
            noalias output: []T,
        ) Error!usize {
            const buf = ctx.typed(T);

            return switch (T) {
                i128, u128 => Impl128(T).bitpack_decompress(
                    buf[0..1024],
                    buf[1024..1536],
                    input,
                    output,
                ),
                i256, u256 => Impl256(T).bitpack_decompress(
                    buf[0..1024],
                    buf[1024..1280],
                    input,
                    output,
                ),
                else => Impl(T).bitpack_decompress(
                    buf[0..1024],
                    input,
                    output,
                ),
            };
        }

        /// Compression function will expect this amount of bytes of capacity on the
        /// output buffer.
        pub fn forpack_compress_bound(len: u32) u32 {
            return switch (T) {
                i128, u128 => Impl128(T).forpack_compress_bound(len),
                i256, u256 => Impl256(T).forpack_compress_bound(len),
                else => Impl(T).forpack_compress_bound(len),
            };
        }

        /// Compress integers from `input` to `output` buffer.
        ///
        /// Will error if `output.len < forpack_compress_bound(input.len)` or if
        /// `input.len >= maxInt(u32)`.
        ///
        /// Returns the number of bytes written to the `output`.
        pub fn forpack_compress(
            ctx: Ctx,
            noalias input: []const T,
            noalias output: []u8,
        ) Error!usize {
            const buf = ctx.typed(T);

            return switch (T) {
                i128, u128 => Impl128(T).forpack_compress(
                    buf[0..1024],
                    buf[1024..1536],
                    input,
                    output,
                ),
                i256, u256 => Impl256(T).forpack_compress(
                    buf[0..1024],
                    buf[1024..1280],
                    input,
                    output,
                ),
                else => Impl(T).forpack_compress(
                    buf[0..1024],
                    input,
                    output,
                ),
            };
        }

        /// Decompress integers from `input` buffer to `output`.
        ///
        /// Expects enough data to decompress `output.len` integers from `input`.
        ///
        /// Returns the byte count of data read from `input` buffer.
        ///
        /// Will error if input data is invalid or short.
        pub fn forpack_decompress(
            ctx: Ctx,
            noalias input: []const u8,
            noalias output: []T,
        ) Error!usize {
            const buf = ctx.typed(T);

            return switch (T) {
                i128, u128 => Impl128(T).forpack_decompress(
                    buf[0..1024],
                    buf[1024..1536],
                    input,
                    output,
                ),
                i256, u256 => Impl256(T).forpack_decompress(
                    buf[0..1024],
                    buf[1024..1280],
                    input,
                    output,
                ),
                else => Impl(T).forpack_decompress(
                    buf[0..1024],
                    input,
                    output,
                ),
            };
        }

        /// Compression function will expect this amount of bytes of capacity on the
        /// output buffer.
        pub fn deltapack_compress_bound(len: u32) u32 {
            return switch (T) {
                i128, u128 => Impl128(T).delta_compress_bound(len),
                i256, u256 => Impl256(T).delta_compress_bound(len),
                else => Impl(T).delta_compress_bound(len),
            };
        }

        /// Compress integers from `input` to `output` buffer.
        ///
        /// Will error if `output.len < deltapack_compress_bound(input.len)` or if
        /// `input.len >= maxInt(u32)`.
        ///
        /// Returns the number of bytes written to the `output`.
        pub fn deltapack_compress(
            ctx: Ctx,
            noalias input: []const T,
            noalias output: []u8,
        ) Error!usize {
            const buf = ctx.typed(T);

            return switch (T) {
                i128, u128 => Impl128(T).delta_compress(
                    buf[0..1024],
                    buf[1024..1536],
                    buf[1536..2048],
                    input,
                    output,
                ),
                i256, u256 => Impl256(T).delta_compress(
                    buf[0..1024],
                    buf[1024..1280],
                    buf[1280..1536],
                    input,
                    output,
                ),
                else => Impl(T).delta_compress(
                    buf[0..1024],
                    buf[1024..2048],
                    input,
                    output,
                ),
            };
        }

        /// Decompress integers from `input` buffer to `output`.
        ///
        /// Expects enough data to decompress `output.len` integers from `input`.
        ///
        /// Returns the byte count of data read from `input` buffer.
        ///
        /// Will error if input data is invalid or short.
        pub fn deltapack_decompress(
            ctx: Ctx,
            noalias input: []const u8,
            noalias output: []T,
        ) Error!usize {
            const buf = ctx.typed(T);

            return switch (T) {
                i128, u128 => Impl128(T).delta_decompress(
                    buf[0..1024],
                    buf[1024..1536],
                    buf[1536..2048],
                    input,
                    output,
                ),
                i256, u256 => Impl256(T).delta_decompress(
                    buf[0..1024],
                    buf[1024..1280],
                    buf[1280..1536],
                    input,
                    output,
                ),
                else => Impl(T).delta_decompress(
                    buf[0..1024],
                    buf[1024..2048],
                    input,
                    output,
                ),
            };
        }
    };
}

fn Impl256(comptime T: type) type {
    const I = switch (T) {
        i256 => i64,
        u256 => u64,
        else => @compileError("unexpected type"),
    };

    const Inner = Impl(I);

    return struct {
        fn bitpack_compress_bound(len: u32) u32 {
            const n_blocks = len / 1024;
            const n_remainder = len % 1024;

            // We will compress the four blocks together as they won't effect each other.
            const size_per_block = Inner.bitpack_compress_bound(4096);

            // We will make four seperate calls to compress the remainder
            const size_for_remainder = 4 * Inner.bitpack_compress_bound(n_remainder);

            return n_blocks * size_per_block + size_for_remainder;
        }

        fn bitpack_compress(
            noalias split_buf: *[1024]T,
            noalias scratch: *[256]T,
            noalias input: []const T,
            noalias output: []u8,
        ) Error!usize {
            if (input.len > std.math.maxInt(u32)) {
                return Error.InvalidInput;
            }
            const len: u32 = @intCast(input.len);

            const n_remainder = len % 1024;

            const output_bound = bitpack_compress_bound(len);
            if (output.len < output_bound) {
                return Error.InvalidInput;
            }

            // number of bytes written so far
            var offset: usize = 0;

            const split_a: *[1024]I = @ptrCast(split_buf[0..256]);
            const split_b: *[1024]I = @ptrCast(split_buf[256..512]);
            const split_c: *[1024]I = @ptrCast(split_buf[512..768]);
            const split_d: *[1024]I = @ptrCast(split_buf[768..1024]);

            const scratch_buf: *[1024]I = @ptrCast(scratch);

            // write remainder data
            if (n_remainder > 0) {
                split(
                    input[0..n_remainder],
                    split_a[0..n_remainder],
                    split_b[0..n_remainder],
                    split_c[0..n_remainder],
                    split_d[0..n_remainder],
                );

                offset += try Inner.bitpack_compress(scratch_buf, split_a[0..n_remainder], output[offset..]);
                offset += try Inner.bitpack_compress(scratch_buf, split_b[0..n_remainder], output[offset..]);
                offset += try Inner.bitpack_compress(scratch_buf, split_c[0..n_remainder], output[offset..]);
                offset += try Inner.bitpack_compress(scratch_buf, split_d[0..n_remainder], output[offset..]);
            }

            const blocks: []const [1024]T = @ptrCast(input[n_remainder..]);
            for (blocks) |*block| {
                split(block, split_a, split_b, split_c, split_d);
                offset += try Inner.bitpack_compress(scratch_buf, @ptrCast(split_buf), output[offset..]);
            }

            return offset;
        }

        fn bitpack_decompress(
            noalias split_buf: *[1024]T,
            noalias scratch: *[256]T,
            noalias input: []const u8,
            noalias output: []T,
        ) Error!usize {
            if (output.len > std.math.maxInt(u32)) {
                return Error.InvalidInput;
            }
            const len: u32 = @intCast(output.len);

            const n_whole_blocks = len / 1024;
            const n_remainder = len % 1024;

            // number of bytes read so far
            var offset: usize = 0;

            const split_a: *[1024]I = @ptrCast(split_buf[0..256]);
            const split_b: *[1024]I = @ptrCast(split_buf[256..512]);
            const split_c: *[1024]I = @ptrCast(split_buf[512..768]);
            const split_d: *[1024]I = @ptrCast(split_buf[768..1024]);

            const scratch_buf: *[1024]I = @ptrCast(scratch);

            if (n_remainder > 0) {
                offset += try Inner.bitpack_decompress(scratch_buf, input[offset..], split_a[0..n_remainder]);
                offset += try Inner.bitpack_decompress(scratch_buf, input[offset..], split_b[0..n_remainder]);
                offset += try Inner.bitpack_decompress(scratch_buf, input[offset..], split_c[0..n_remainder]);
                offset += try Inner.bitpack_decompress(scratch_buf, input[offset..], split_d[0..n_remainder]);

                combine(
                    split_a[0..n_remainder],
                    split_b[0..n_remainder],
                    split_c[0..n_remainder],
                    split_d[0..n_remainder],
                    output[0..n_remainder],
                );
            }

            for (0..n_whole_blocks) |block_idx| {
                offset += try Inner.bitpack_decompress(scratch_buf, input[offset..], @ptrCast(split_buf));

                const out_offset = output[block_idx * 1024 + n_remainder ..];
                combine(split_a, split_b, split_c, split_d, out_offset[0..1024]);
            }

            return offset;
        }

        fn forpack_compress_bound(len: u32) u32 {
            const n_blocks = len / 1024;
            const n_remainder = len % 1024;

            // We will compress the four blocks together as they won't effect each other.
            const size_per_block = Inner.forpack_compress_bound(4096);

            // We will make four seperate calls to compress the remainder
            const size_for_remainder = 4 * Inner.forpack_compress_bound(n_remainder);

            return n_blocks * size_per_block + size_for_remainder;
        }

        fn forpack_compress(
            noalias split_buf: *[1024]T,
            noalias scratch: *[256]T,
            noalias input: []const T,
            noalias output: []u8,
        ) Error!usize {
            if (input.len > std.math.maxInt(u32)) {
                return Error.InvalidInput;
            }
            const len: u32 = @intCast(input.len);

            const n_remainder = len % 1024;

            const output_bound = forpack_compress_bound(len);
            if (output.len < output_bound) {
                return Error.InvalidInput;
            }

            // number of bytes written so far
            var offset: usize = 0;

            const split_a: *[1024]I = @ptrCast(split_buf[0..256]);
            const split_b: *[1024]I = @ptrCast(split_buf[256..512]);
            const split_c: *[1024]I = @ptrCast(split_buf[512..768]);
            const split_d: *[1024]I = @ptrCast(split_buf[768..1024]);

            const scratch_buf: *[1024]I = @ptrCast(scratch);

            // write remainder data
            if (n_remainder > 0) {
                split(
                    input[0..n_remainder],
                    split_a[0..n_remainder],
                    split_b[0..n_remainder],
                    split_c[0..n_remainder],
                    split_d[0..n_remainder],
                );

                offset += try Inner.forpack_compress(scratch_buf, split_a[0..n_remainder], output[offset..]);
                offset += try Inner.forpack_compress(scratch_buf, split_b[0..n_remainder], output[offset..]);
                offset += try Inner.forpack_compress(scratch_buf, split_c[0..n_remainder], output[offset..]);
                offset += try Inner.forpack_compress(scratch_buf, split_d[0..n_remainder], output[offset..]);
            }

            const blocks: []const [1024]T = @ptrCast(input[n_remainder..]);
            for (blocks) |*block| {
                split(block, split_a, split_b, split_c, split_d);
                offset += try Inner.forpack_compress(scratch_buf, @ptrCast(split_buf), output[offset..]);
            }

            return offset;
        }

        fn forpack_decompress(
            noalias split_buf: *[1024]T,
            noalias scratch: *[256]T,
            noalias input: []const u8,
            noalias output: []T,
        ) Error!usize {
            if (output.len > std.math.maxInt(u32)) {
                return Error.InvalidInput;
            }
            const len: u32 = @intCast(output.len);

            const n_whole_blocks = len / 1024;
            const n_remainder = len % 1024;

            // number of bytes read so far
            var offset: usize = 0;

            const split_a: *[1024]I = @ptrCast(split_buf[0..256]);
            const split_b: *[1024]I = @ptrCast(split_buf[256..512]);
            const split_c: *[1024]I = @ptrCast(split_buf[512..768]);
            const split_d: *[1024]I = @ptrCast(split_buf[768..1024]);

            const scratch_buf: *[1024]I = @ptrCast(scratch);

            if (n_remainder > 0) {
                offset += try Inner.forpack_decompress(scratch_buf, input[offset..], split_a[0..n_remainder]);
                offset += try Inner.forpack_decompress(scratch_buf, input[offset..], split_b[0..n_remainder]);
                offset += try Inner.forpack_decompress(scratch_buf, input[offset..], split_c[0..n_remainder]);
                offset += try Inner.forpack_decompress(scratch_buf, input[offset..], split_d[0..n_remainder]);

                combine(
                    split_a[0..n_remainder],
                    split_b[0..n_remainder],
                    split_c[0..n_remainder],
                    split_d[0..n_remainder],
                    output[0..n_remainder],
                );
            }

            for (0..n_whole_blocks) |block_idx| {
                offset += try Inner.forpack_decompress(scratch_buf, input[offset..], @ptrCast(split_buf));

                const out_offset = output[block_idx * 1024 + n_remainder ..];
                combine(
                    split_a,
                    split_b,
                    split_c,
                    split_d,
                    out_offset[0..1024],
                );
            }

            return offset;
        }

        fn delta_compress_bound(len: u32) u32 {
            const n_blocks = len / 1024;
            const n_remainder = len % 1024;

            // We will compress the four blocks together as they won't effect each other.
            const size_per_block = Inner.delta_compress_bound(4096);

            // We will make four seperate calls to compress the remainder
            const size_for_remainder = 4 * Inner.delta_compress_bound(n_remainder);

            return n_blocks * size_per_block + size_for_remainder;
        }

        fn delta_compress(
            noalias split_buf: *[1024]T,
            noalias transposed: *[256]T,
            noalias delta: *[256]T,
            noalias input: []const T,
            noalias output: []u8,
        ) Error!usize {
            if (input.len > std.math.maxInt(u32)) {
                return Error.InvalidInput;
            }
            const len: u32 = @intCast(input.len);

            const n_remainder = len % 1024;

            const output_bound = delta_compress_bound(len);
            if (output.len < output_bound) {
                return Error.InvalidInput;
            }

            // number of bytes written so far
            var offset: usize = 0;

            const split_a: *[1024]I = @ptrCast(split_buf[0..256]);
            const split_b: *[1024]I = @ptrCast(split_buf[256..512]);
            const split_c: *[1024]I = @ptrCast(split_buf[512..768]);
            const split_d: *[1024]I = @ptrCast(split_buf[768..1024]);

            const transposed_buf: *[1024]I = @ptrCast(transposed);
            const delta_buf: *[1024]I = @ptrCast(delta);

            // write remainder data
            if (n_remainder > 0) {
                split(
                    input[0..n_remainder],
                    split_a[0..n_remainder],
                    split_b[0..n_remainder],
                    split_c[0..n_remainder],
                    split_d[0..n_remainder],
                );

                offset += try Inner.delta_compress(transposed_buf, delta_buf, split_a[0..n_remainder], output[offset..]);
                offset += try Inner.delta_compress(transposed_buf, delta_buf, split_b[0..n_remainder], output[offset..]);
                offset += try Inner.delta_compress(transposed_buf, delta_buf, split_c[0..n_remainder], output[offset..]);
                offset += try Inner.delta_compress(transposed_buf, delta_buf, split_d[0..n_remainder], output[offset..]);
            }

            const blocks: []const [1024]T = @ptrCast(input[n_remainder..]);
            for (blocks) |*block| {
                split(block, split_a, split_b, split_c, split_d);
                offset += try Inner.delta_compress(transposed_buf, delta_buf, @ptrCast(split_buf), output[offset..]);
            }

            return offset;
        }

        fn delta_decompress(
            noalias split_buf: *[1024]T,
            noalias scratch: *[256]T,
            noalias transposed: *[256]T,
            noalias input: []const u8,
            noalias output: []T,
        ) Error!usize {
            if (output.len > std.math.maxInt(u32)) {
                return Error.InvalidInput;
            }
            const len: u32 = @intCast(output.len);

            const n_whole_blocks = len / 1024;
            const n_remainder = len % 1024;

            // number of bytes read so far
            var offset: usize = 0;

            const split_a: *[1024]I = @ptrCast(split_buf[0..256]);
            const split_b: *[1024]I = @ptrCast(split_buf[256..512]);
            const split_c: *[1024]I = @ptrCast(split_buf[512..768]);
            const split_d: *[1024]I = @ptrCast(split_buf[768..1024]);

            const scratch_buf: *[1024]I = @ptrCast(scratch);
            const transposed_buf: *[1024]I = @ptrCast(transposed);

            if (n_remainder > 0) {
                offset += try Inner.delta_decompress(
                    scratch_buf,
                    transposed_buf,
                    input[offset..],
                    split_a[0..n_remainder],
                );
                offset += try Inner.delta_decompress(
                    scratch_buf,
                    transposed_buf,
                    input[offset..],
                    split_b[0..n_remainder],
                );
                offset += try Inner.delta_decompress(
                    scratch_buf,
                    transposed_buf,
                    input[offset..],
                    split_c[0..n_remainder],
                );
                offset += try Inner.delta_decompress(
                    scratch_buf,
                    transposed_buf,
                    input[offset..],
                    split_d[0..n_remainder],
                );

                combine(
                    split_a[0..n_remainder],
                    split_b[0..n_remainder],
                    split_c[0..n_remainder],
                    split_d[0..n_remainder],
                    output[0..n_remainder],
                );
            }

            for (0..n_whole_blocks) |block_idx| {
                offset += try Inner.delta_decompress(
                    scratch_buf,
                    transposed_buf,
                    input[offset..],
                    @ptrCast(split_buf),
                );

                const out_offset = output[block_idx * 1024 + n_remainder ..];
                combine(
                    split_a,
                    split_b,
                    split_c,
                    split_d,
                    out_offset[0..1024],
                );
            }

            return offset;
        }

        fn combine1024(
            noalias a: *const [1024]I,
            noalias b: *const [1024]I,
            noalias c: *const [1024]I,
            noalias d: *const [1024]I,
            noalias output: *[1024]T,
        ) void {
            const out: *[1024][4]I = @ptrCast(output);
            for (0..1024) |i| {
                out[i] = .{ a[i], b[i], c[i], d[i] };
            }
        }

        fn combine(
            noalias a: []const I,
            noalias b: []const I,
            noalias c: []const I,
            noalias d: []const I,
            noalias output: []T,
        ) void {
            std.debug.assert(output.len == a.len);
            std.debug.assert(output.len == b.len);
            std.debug.assert(output.len == c.len);
            std.debug.assert(output.len == d.len);

            const out: [][4]I = @ptrCast(output);
            for (out, a, b, c, d) |*o, a_, b_, c_, d_| {
                o.* = .{ a_, b_, c_, d_ };
            }
        }

        fn split1024(
            noalias input: *const [1024]T,
            noalias a: *[1024]I,
            noalias b: *[1024]I,
            noalias c: *[1024]I,
            noalias d: *[1024]I,
        ) void {
            for (0..1024) |i| {
                const v: [4]I = @bitCast(input[i]);
                a[i] = v[0];
                b[i] = v[1];
                c[i] = v[2];
                d[i] = v[3];
            }
        }

        fn split(
            noalias input: []const T,
            noalias a: []I,
            noalias b: []I,
            noalias c: []I,
            noalias d: []I,
        ) void {
            std.debug.assert(input.len == a.len);
            std.debug.assert(input.len == b.len);
            std.debug.assert(input.len == c.len);
            std.debug.assert(input.len == d.len);

            for (input, a, b, c, d) |i, *a_, *b_, *c_, *d_| {
                const v: [4]I = @bitCast(i);
                a_.* = v[0];
                b_.* = v[1];
                c_.* = v[2];
                d_.* = v[3];
            }
        }
    };
}

fn Impl128(comptime T: type) type {
    const I = switch (T) {
        i128 => i64,
        u128 => u64,
        else => @compileError("unexpected type"),
    };

    const Inner = Impl(I);

    return struct {
        fn bitpack_compress_bound(len: u32) u32 {
            const n_blocks = len / 1024;
            const n_remainder = len % 1024;

            // We will compress the two blocks together as they won't effect each other.
            const size_per_block = Inner.bitpack_compress_bound(2048);

            // We will make two seperate calls to compress the remainder
            const size_for_remainder = 2 * Inner.bitpack_compress_bound(n_remainder);

            return n_blocks * size_per_block + size_for_remainder;
        }

        fn bitpack_compress(
            noalias split_buf: *[1024]T,
            noalias scratch: *[512]T,
            noalias input: []const T,
            noalias output: []u8,
        ) Error!usize {
            if (input.len > std.math.maxInt(u32)) {
                return Error.InvalidInput;
            }
            const len: u32 = @intCast(input.len);

            const n_remainder = len % 1024;

            const output_bound = bitpack_compress_bound(len);
            if (output.len < output_bound) {
                return Error.InvalidInput;
            }

            // number of bytes written so far
            var offset: usize = 0;

            const split_a: *[1024]I = @ptrCast(split_buf[0..512]);
            const split_b: *[1024]I = @ptrCast(split_buf[512..]);
            const scratch_buf: *[1024]I = @ptrCast(scratch);

            // write remainder data
            if (n_remainder > 0) {
                split(input[0..n_remainder], split_a[0..n_remainder], split_b[0..n_remainder]);

                offset += try Inner.bitpack_compress(scratch_buf, split_a[0..n_remainder], output[offset..]);
                offset += try Inner.bitpack_compress(scratch_buf, split_b[0..n_remainder], output[offset..]);
            }

            const blocks: []const [1024]T = @ptrCast(input[n_remainder..]);
            for (blocks) |*block| {
                split(block, split_a, split_b);
                offset += try Inner.bitpack_compress(scratch_buf, @ptrCast(split_buf), output[offset..]);
            }

            return offset;
        }

        fn bitpack_decompress(
            noalias split_buf: *[1024]T,
            noalias scratch: *[512]T,
            noalias input: []const u8,
            noalias output: []T,
        ) Error!usize {
            if (output.len > std.math.maxInt(u32)) {
                return Error.InvalidInput;
            }
            const len: u32 = @intCast(output.len);

            const n_whole_blocks = len / 1024;
            const n_remainder = len % 1024;

            // number of bytes read so far
            var offset: usize = 0;

            const split_a: *[1024]I = @ptrCast(split_buf[0..512]);
            const split_b: *[1024]I = @ptrCast(split_buf[512..]);
            const scratch_buf: *[1024]I = @ptrCast(scratch);

            if (n_remainder > 0) {
                offset += try Inner.bitpack_decompress(scratch_buf, input[offset..], split_a[0..n_remainder]);
                offset += try Inner.bitpack_decompress(scratch_buf, input[offset..], split_b[0..n_remainder]);

                combine(split_a[0..n_remainder], split_b[0..n_remainder], output[0..n_remainder]);
            }

            for (0..n_whole_blocks) |block_idx| {
                offset += try Inner.bitpack_decompress(scratch_buf, input[offset..], @ptrCast(split_buf));

                const out_offset = output[block_idx * 1024 + n_remainder ..];
                combine(split_a, split_b, out_offset[0..1024]);
            }

            return offset;
        }

        fn forpack_compress_bound(len: u32) u32 {
            const n_blocks = len / 1024;
            const n_remainder = len % 1024;

            // We will compress the two blocks together as they won't effect each other.
            const size_per_block = Inner.forpack_compress_bound(2048);

            // We will make two seperate calls to compress the remainder
            const size_for_remainder = 2 * Inner.forpack_compress_bound(n_remainder);

            return n_blocks * size_per_block + size_for_remainder;
        }

        fn forpack_compress(
            noalias split_buf: *[1024]T,
            noalias scratch: *[512]T,
            noalias input: []const T,
            noalias output: []u8,
        ) Error!usize {
            if (input.len > std.math.maxInt(u32)) {
                return Error.InvalidInput;
            }
            const len: u32 = @intCast(input.len);

            const n_remainder = len % 1024;

            const output_bound = forpack_compress_bound(len);
            if (output.len < output_bound) {
                return Error.InvalidInput;
            }

            // number of bytes written so far
            var offset: usize = 0;

            const split_a: *[1024]I = @ptrCast(split_buf[0..512]);
            const split_b: *[1024]I = @ptrCast(split_buf[512..]);
            const scratch_buf: *[1024]I = @ptrCast(scratch);

            // write remainder data
            if (n_remainder > 0) {
                split(input[0..n_remainder], split_a[0..n_remainder], split_b[0..n_remainder]);

                offset += try Inner.forpack_compress(scratch_buf, split_a[0..n_remainder], output[offset..]);
                offset += try Inner.forpack_compress(scratch_buf, split_b[0..n_remainder], output[offset..]);
            }

            const blocks: []const [1024]T = @ptrCast(input[n_remainder..]);
            for (blocks) |*block| {
                split(block, split_a, split_b);
                offset += try Inner.forpack_compress(scratch_buf, @ptrCast(split_buf), output[offset..]);
            }

            return offset;
        }

        fn forpack_decompress(
            noalias split_buf: *[1024]T,
            noalias scratch: *[512]T,
            noalias input: []const u8,
            noalias output: []T,
        ) Error!usize {
            if (output.len > std.math.maxInt(u32)) {
                return Error.InvalidInput;
            }
            const len: u32 = @intCast(output.len);

            const n_whole_blocks = len / 1024;
            const n_remainder = len % 1024;

            // number of bytes read so far
            var offset: usize = 0;

            const split_a: *[1024]I = @ptrCast(split_buf[0..512]);
            const split_b: *[1024]I = @ptrCast(split_buf[512..]);
            const scratch_buf: *[1024]I = @ptrCast(scratch);

            if (n_remainder > 0) {
                offset += try Inner.forpack_decompress(scratch_buf, input[offset..], split_a[0..n_remainder]);
                offset += try Inner.forpack_decompress(scratch_buf, input[offset..], split_b[0..n_remainder]);

                combine(split_a[0..n_remainder], split_b[0..n_remainder], output[0..n_remainder]);
            }

            for (0..n_whole_blocks) |block_idx| {
                offset += try Inner.forpack_decompress(scratch_buf, input[offset..], @ptrCast(split_buf));

                const out_offset = output[block_idx * 1024 + n_remainder ..];
                combine(split_a, split_b, out_offset[0..1024]);
            }

            return offset;
        }

        fn delta_compress_bound(len: u32) u32 {
            const n_blocks = len / 1024;
            const n_remainder = len % 1024;

            // We will compress the two blocks together as they won't effect each other.
            const size_per_block = Inner.delta_compress_bound(2048);

            // We will make two seperate calls to compress the remainder
            const size_for_remainder = 2 * Inner.delta_compress_bound(n_remainder);

            return n_blocks * size_per_block + size_for_remainder;
        }

        fn delta_compress(
            noalias split_buf: *[1024]T,
            noalias transposed: *[512]T,
            noalias delta: *[512]T,
            noalias input: []const T,
            noalias output: []u8,
        ) Error!usize {
            if (input.len > std.math.maxInt(u32)) {
                return Error.InvalidInput;
            }
            const len: u32 = @intCast(input.len);

            const n_remainder = len % 1024;

            const output_bound = delta_compress_bound(len);
            if (output.len < output_bound) {
                return Error.InvalidInput;
            }

            // number of bytes written so far
            var offset: usize = 0;

            const split_a: *[1024]I = @ptrCast(split_buf[0..512]);
            const split_b: *[1024]I = @ptrCast(split_buf[512..]);
            const transposed_buf: *[1024]I = @ptrCast(transposed);
            const delta_buf: *[1024]I = @ptrCast(delta);

            // write remainder data
            if (n_remainder > 0) {
                split(input[0..n_remainder], split_a[0..n_remainder], split_b[0..n_remainder]);

                offset += try Inner.delta_compress(transposed_buf, delta_buf, split_a[0..n_remainder], output[offset..]);
                offset += try Inner.delta_compress(transposed_buf, delta_buf, split_b[0..n_remainder], output[offset..]);
            }

            const blocks: []const [1024]T = @ptrCast(input[n_remainder..]);
            for (blocks) |*block| {
                split(block, split_a, split_b);
                offset += try Inner.delta_compress(transposed_buf, delta_buf, @ptrCast(split_buf), output[offset..]);
            }

            return offset;
        }

        fn delta_decompress(
            noalias split_buf: *[1024]T,
            noalias scratch: *[512]T,
            noalias transposed: *[512]T,
            noalias input: []const u8,
            noalias output: []T,
        ) Error!usize {
            if (output.len > std.math.maxInt(u32)) {
                return Error.InvalidInput;
            }
            const len: u32 = @intCast(output.len);

            const n_whole_blocks = len / 1024;
            const n_remainder = len % 1024;

            // number of bytes read so far
            var offset: usize = 0;

            const split_a: *[1024]I = @ptrCast(split_buf[0..512]);
            const split_b: *[1024]I = @ptrCast(split_buf[512..]);
            const scratch_buf: *[1024]I = @ptrCast(scratch);
            const transposed_buf: *[1024]I = @ptrCast(transposed);

            if (n_remainder > 0) {
                offset += try Inner.delta_decompress(
                    scratch_buf,
                    transposed_buf,
                    input[offset..],
                    split_a[0..n_remainder],
                );
                offset += try Inner.delta_decompress(
                    scratch_buf,
                    transposed_buf,
                    input[offset..],
                    split_b[0..n_remainder],
                );

                combine(split_a[0..n_remainder], split_b[0..n_remainder], output[0..n_remainder]);
            }

            for (0..n_whole_blocks) |block_idx| {
                offset += try Inner.delta_decompress(
                    scratch_buf,
                    transposed_buf,
                    input[offset..],
                    @ptrCast(split_buf),
                );

                const out_offset = output[block_idx * 1024 + n_remainder ..];
                combine(split_a, split_b, out_offset[0..1024]);
            }

            return offset;
        }

        fn combine1024(
            noalias a: *const [1024]I,
            noalias b: *const [1024]I,
            noalias output: *[1024]T,
        ) void {
            const out: *[1024][2]I = @ptrCast(output);
            for (0..1024) |i| {
                out[i] = .{ a[i], b[i] };
            }
        }

        fn combine(noalias a: []const I, noalias b: []const I, noalias output: []T) void {
            std.debug.assert(output.len == a.len);
            std.debug.assert(output.len == b.len);

            const out: [][2]I = @ptrCast(output);
            for (out, a, b) |*o, a_, b_| {
                o.* = .{ a_, b_ };
            }
        }

        fn split1024(noalias input: *const [1024]T, noalias a: *[1024]I, noalias b: *[1024]I) void {
            for (0..1024) |i| {
                const v: [2]I = @bitCast(input[i]);
                a[i] = v[0];
                b[i] = v[1];
            }
        }

        fn split(noalias input: []const T, noalias a: []I, noalias b: []I) void {
            std.debug.assert(input.len == a.len);
            std.debug.assert(input.len == b.len);

            for (input, a, b) |i, *a_, *b_| {
                const v: [2]I = @bitCast(i);
                a_.* = v[0];
                b_.* = v[1];
            }
        }
    };
}

fn Impl(comptime T: type) type {
    const U = switch (T) {
        u8, u16, u32, u64 => T,
        i8 => u8,
        i16 => u16,
        i32 => u32,
        i64 => u64,
        else => @compileError("unexpected type."),
    };

    const IS_SIGNED = T != U;

    const FL = FastLanes(U);
    const SPack = ScalarBitpack(U);

    const N_BYTES = @sizeOf(U);
    const N_BITS = N_BYTES * 8;

    return struct {
        /// Returns the number of BYTES needed on the output buffer for compressing
        /// `len` elements. The compression will likely end up writing less data
        /// but this is needed to make sure the compression will succeed.
        fn bitpack_compress_bound(len: u32) u32 {
            const n_blocks = len / 1024;
            const n_remainder = len % 1024;

            const max_block_size = N_BYTES * 1024;

            // Layout of the output is:
            // - remainder values bit_width(u8)
            // - bit_width(u8) per block
            // - packed remainder values
            // - packed full blocks

            return 1 + n_blocks + n_remainder * N_BYTES + max_block_size * n_blocks;
        }

        /// Compress the input integers into the output buffer.
        /// `output` should be at least `bitpack_compress_bound(input.len)` BYTES.
        ///
        /// Returns the number of bytes written to the output.
        fn bitpack_compress(
            noalias scratch: *[1024]T,
            noalias input: []const T,
            noalias output: []u8,
        ) Error!usize {
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

            // Start writing compressed data, skip some bytes for saving bit_widths later on
            // 1 is for the byte_width of remainder data, rest is for
            // byte widths of whole blocks.
            const byte_offset = 1 + n_whole_blocks;

            // We should have this much capacity even though we will write less
            const out_t_len = 1024 * n_whole_blocks + n_remainder;

            const out: []align(1) U = @ptrCast(output[byte_offset .. byte_offset + out_t_len * N_BYTES]);
            std.debug.assert(out.len == out_t_len);

            // number of T written so far, NOT number of bytes
            var offset: usize = 0;

            // write remainder data
            if (n_remainder > 0) {
                const remainder_data = load_remainder(
                    input[0..n_remainder],
                    scratch[0..n_remainder],
                );

                const max = std.mem.max(U, remainder_data);
                const width = needed_width(max);

                output[0] = width;

                const remainder_packed_len = (@as(u64, width) * n_remainder + N_BITS - 1) / N_BITS;
                const n_written = SPack.bitpack(
                    remainder_data,
                    0,
                    out,
                    width,
                ) catch {
                    @panic("failed scalar pack, this should never happen");
                };
                std.debug.assert(n_written == remainder_packed_len);
                offset += remainder_packed_len;
            } else {
                output[0] = 0;
            }

            // Write whole blocks
            const input_blocks: []const [1024]T = @ptrCast(input[n_remainder..]);
            for (0..n_whole_blocks) |block_idx| {
                const block = load_block(&input_blocks[block_idx], scratch);

                const max = max1024(block);
                const width = needed_width(max);

                output[1 + block_idx] = width;

                offset += FL.dyn_bit_pack(block, out[offset..], width);
            }

            return byte_offset + N_BYTES * offset;
        }

        /// Returns the number of bytes decompressed from the `input`
        fn bitpack_decompress(
            noalias scratch: *[1024]T,
            noalias input: []const u8,
            noalias output: []T,
        ) Error!usize {
            if (output.len > std.math.maxInt(u32)) {
                return Error.InvalidInput;
            }
            const len: u32 = @intCast(output.len);

            const n_whole_blocks = len / 1024;
            const n_remainder = len % 1024;

            if (input.len < 1 + n_whole_blocks) {
                return Error.InvalidInput;
            }

            const remainder_width = input[0];
            if (remainder_width > N_BITS) {
                return Error.InvalidInput;
            }

            const block_widths = input[1 .. 1 + n_whole_blocks];

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

            const data_section_byte_offset = 1 + n_whole_blocks;

            const data_section_bytes = input[data_section_byte_offset..];
            if (data_section_bytes.len < total_packed_len * N_BYTES) {
                return Error.InvalidInput;
            }

            const data_section: []align(1) const U = @ptrCast(data_section_bytes[0 .. total_packed_len * N_BYTES]);

            // number of T read so far, NOT number of bytes
            var offset: usize = 0;

            // read remainder data
            if (n_remainder > 0) {
                if (IS_SIGNED) {
                    const zigzagged: []U = @ptrCast(scratch[0..n_remainder]);
                    const n_read = try SPack.bitunpack(
                        data_section[0..remainder_packed_len],
                        0,
                        zigzagged,
                        remainder_width,
                    );
                    std.debug.assert(n_read == remainder_packed_len);
                    ZigZag(T).decode(zigzagged, output[0..n_remainder]);
                } else {
                    const scratch_u: []U = @ptrCast(scratch[0..n_remainder]);
                    const n_read = try SPack.bitunpack(
                        data_section[0..remainder_packed_len],
                        0,
                        scratch_u,
                        remainder_width,
                    );
                    std.debug.assert(n_read == remainder_packed_len);
                    @memcpy(output[0..n_remainder], scratch_u[0..n_remainder]);
                }
                offset += remainder_packed_len;
            }

            // Read whole blocks
            for (0..n_whole_blocks) |block_idx| {
                const width = block_widths[block_idx];

                const out_offset = output[block_idx * 1024 + n_remainder ..];
                if (IS_SIGNED) {
                    const zigzagged: *[1024]U = @ptrCast(scratch);
                    offset += FL.dyn_bit_unpack(
                        data_section[offset..],
                        zigzagged,
                        width,
                    );
                    ZigZag(T).decode1024(zigzagged, out_offset[0..1024]);
                } else {
                    const scratch_u: *[1024]U = @ptrCast(scratch);
                    offset += FL.dyn_bit_unpack(
                        data_section[offset..],
                        scratch_u,
                        width,
                    );
                    out_offset[0..1024].* = scratch_u.*;
                }
            }

            std.debug.assert(offset == total_packed_len);
            std.debug.assert(offset == data_section.len);

            return data_section_byte_offset + @sizeOf(T) * offset;
        }

        /// Returns the number of BYTES needed on the output buffer for compressing
        /// `len` elements. The compression will likely end up writing less data
        /// but this is needed to make sure the compression will succeed.
        fn forpack_compress_bound(len: u32) u32 {
            const n_blocks = len / 1024;
            const n_remainder = len % 1024;

            const max_block_size = N_BYTES * (1024 + 1);

            // Layout of the output is:
            // - remainder values bit_width(u8)
            // - bit_width(u8) per block
            // - packed remainder values, a reference value before the packed values
            // - packed full blocks, a reference value before every block

            return 1 + n_blocks + (1 + n_remainder) * N_BYTES + max_block_size * n_blocks;
        }

        /// Compress the input integers into the output buffer.
        /// `output` should be at least `forpack_compress_bound(input.len)` BYTES.
        ///
        /// Returns the number of bytes written to the output.
        fn forpack_compress(
            noalias scratch: *[1024]T,
            noalias input: []const T,
            noalias output: []u8,
        ) Error!usize {
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

            // Start writing compressed data, skip some bytes for saving bit_widths later on
            // 1 is for the byte_width of remainder data, rest is for
            // byte widths of whole blocks.
            const byte_offset = 1 + n_whole_blocks;

            // We should have this much capacity even though we will write less
            const out_t_len = (1024 + 1) * n_whole_blocks + 1 + n_remainder;

            const out: []align(1) U = @ptrCast(output[byte_offset .. byte_offset + out_t_len * N_BYTES]);
            std.debug.assert(out.len == out_t_len);

            // number of T written so far, NOT number of bytes
            var offset: usize = 0;

            // write remainder data
            if (n_remainder > 0) {
                const remainder_data = load_remainder(input[0..n_remainder], scratch[0..n_remainder]);

                const min, const max = std.mem.minMax(U, remainder_data);
                const width = needed_width(max - min);

                output[0] = width;

                out[0] = min;
                offset += 1;

                const remainder_packed_len = (@as(u64, width) * n_remainder + N_BITS - 1) / N_BITS;
                const n_written = SPack.bitpack(
                    remainder_data,
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

                output[0] = 0;
            }

            // Write whole blocks
            const input_blocks: []const [1024]T = @ptrCast(input[n_remainder..]);
            for (0..n_whole_blocks) |block_idx| {
                const block = load_block(&input_blocks[block_idx], scratch);

                const min, const max = minmax1024(block);
                const width = needed_width(max - min);

                output[1 + block_idx] = width;

                out[offset] = min;
                offset += 1;

                offset += FL.dyn_for_pack(block, min, out[offset..], width);
            }

            return byte_offset + N_BYTES * offset;
        }

        /// Returns the number of bytes decompressed from the `input`
        fn forpack_decompress(
            noalias scratch: *[1024]T,
            noalias input: []const u8,
            noalias output: []T,
        ) Error!usize {
            if (output.len > std.math.maxInt(u32)) {
                return Error.InvalidInput;
            }
            const len: u32 = @intCast(output.len);

            const n_whole_blocks = len / 1024;
            const n_remainder = len % 1024;

            if (input.len < 1 + n_whole_blocks) {
                return Error.InvalidInput;
            }

            const remainder_width = input[0];
            if (remainder_width > N_BITS) {
                return Error.InvalidInput;
            }

            const block_widths = input[1 .. 1 + n_whole_blocks];

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

            const data_section_byte_offset = 1 + n_whole_blocks;

            const data_section_bytes = input[data_section_byte_offset..];
            if (data_section_bytes.len < total_packed_len * N_BYTES) {
                return Error.InvalidInput;
            }

            const data_section: []align(1) const U = @ptrCast(data_section_bytes[0 .. total_packed_len * N_BYTES]);

            // number of T read so far, NOT number of bytes
            var offset: usize = 0;

            // read remainder data
            if (n_remainder > 0) {
                const ref = data_section[offset];
                offset += 1;

                if (IS_SIGNED) {
                    const zigzagged: []U = @ptrCast(scratch[0..n_remainder]);
                    const n_read = try SPack.bitunpack(
                        data_section[offset .. offset + remainder_packed_len],
                        ref,
                        zigzagged,
                        remainder_width,
                    );
                    std.debug.assert(n_read == remainder_packed_len);
                    ZigZag(T).decode(zigzagged, output[0..n_remainder]);
                } else {
                    const scratch_u: []U = @ptrCast(scratch[0..n_remainder]);
                    const n_read = try SPack.bitunpack(
                        data_section[offset .. offset + remainder_packed_len],
                        ref,
                        scratch_u,
                        remainder_width,
                    );
                    std.debug.assert(n_read == remainder_packed_len);
                    @memcpy(output[0..n_remainder], scratch_u);
                }
                offset += remainder_packed_len;
            } else {
                // as if we read remainder data min value
                offset += 1;
            }

            // Read whole blocks
            for (0..n_whole_blocks) |block_idx| {
                const width = block_widths[block_idx];

                const ref = data_section[offset];
                offset += 1;

                const out_offset = output[block_idx * 1024 + n_remainder ..];
                if (IS_SIGNED) {
                    const zigzagged: *[1024]U = @ptrCast(scratch);
                    offset += FL.dyn_for_unpack(
                        data_section[offset..],
                        ref,
                        zigzagged,
                        width,
                    );
                    ZigZag(T).decode1024(zigzagged, out_offset[0..1024]);
                } else {
                    const scratch_u: *[1024]U = @ptrCast(scratch);
                    offset += FL.dyn_for_unpack(
                        data_section[offset..],
                        ref,
                        scratch_u,
                        width,
                    );
                    out_offset[0..1024].* = scratch_u.*;
                }
            }

            std.debug.assert(offset == total_packed_len);
            std.debug.assert(offset == data_section.len);

            return data_section_byte_offset + @sizeOf(T) * offset;
        }

        /// Returns the number of BYTES needed on the output buffer for compressing
        /// `len` elements. The compression will likely end up writing less data
        /// but this is needed to make sure the compression will succeed.
        fn delta_compress_bound(len: u32) u32 {
            const n_blocks = len / 1024;
            const n_remainder = len % 1024;

            const max_block_size = N_BYTES * (FL.N_LANES + 1024);

            // Layout of the output is:
            // - remainder values bit_width(u8)
            // - bit_width(u8) per block
            // - packed remainder values, a base value before the values
            // - packed full blocks, N_LANES "bases" T before each block

            return 1 + n_blocks + (1 + n_remainder) * N_BYTES + max_block_size * n_blocks;
        }

        /// Compress the input integers into the output buffer.
        /// `output` should be at least `delta_compress_bound(input.len)` BYTES.
        ///
        /// `transposed` and `delta` inputs are for using as internal scratch memory.
        ///
        /// Returns the number of bytes written to the output.
        fn delta_compress(
            noalias transposed: *[1024]T,
            noalias delta: *[1024]T,
            noalias input: []const T,
            noalias output: []u8,
        ) Error!usize {
            if (input.len > std.math.maxInt(u32)) {
                return Error.InvalidInput;
            }
            const len: u32 = @intCast(input.len);

            const n_whole_blocks = len / 1024;
            const n_remainder = len % 1024;

            const output_bound = delta_compress_bound(len);
            if (output.len < output_bound) {
                return Error.InvalidInput;
            }

            // Start writing compressed data, skip some bytes for saving bit_widths later on
            // 1 is for the byte_width of remainder data, rest is for
            // byte widths of whole blocks.
            const byte_offset = 1 + n_whole_blocks;

            // We should have this much capacity even though we will write less
            const out_t_len = (FL.N_LANES + 1024) * n_whole_blocks + 1 + n_remainder;

            const out: []align(1) U = @ptrCast(output[byte_offset .. byte_offset + out_t_len * N_BYTES]);
            std.debug.assert(out.len == out_t_len);

            // number of T written so far, NOT number of bytes
            var offset: usize = 0;

            // write remainder data
            if (n_remainder > 0) {
                // use the delta buffer for zigzag
                const scratch = delta;
                const remainder_data = load_remainder(input[0..n_remainder], scratch[0..n_remainder]);

                const base = remainder_data[0];

                var prev: U = base;
                var max_delta: U = 0;
                for (remainder_data) |v| {
                    const d = v -% prev;
                    max_delta = @max(max_delta, d);
                    prev = v;
                }
                const width = needed_width(max_delta);

                output[0] = width;

                out[offset] = base;
                offset += 1;

                const remainder_packed_len = (@as(u64, width) * n_remainder + N_BITS - 1) / N_BITS;
                const n_written = SPack.delta_bitpack(
                    base,
                    remainder_data,
                    out[offset..],
                    width,
                ) catch {
                    @panic("failed scalar pack, this should never happen");
                };
                std.debug.assert(n_written == remainder_packed_len);
                offset += remainder_packed_len;
            } else {
                output[0] = 0;

                out[offset] = 0;
                offset += 1;
            }

            // Write whole blocks
            const input_blocks: []const [1024]T = @ptrCast(input[n_remainder..]);
            for (0..n_whole_blocks) |block_idx| {
                // use delta for doing zigzag
                // can't use transpose for it because we will write from the zigzag buffer into
                // transpose
                const block = load_block(&input_blocks[block_idx], delta);

                const transposed_buf: *[1024]U = @ptrCast(transposed);
                const delta_buf: *[1024]U = @ptrCast(delta);

                FL.transpose(block, transposed_buf);

                const bases: *const [FL.N_LANES]U = transposed_buf[0..FL.N_LANES];

                FL.delta(transposed_buf, bases, delta_buf);

                const max = max1024(delta_buf);
                const width = needed_width(max);

                output[1 + block_idx] = width;

                for (0..FL.N_LANES) |i| {
                    out[offset + i] = bases[i];
                }
                offset += FL.N_LANES;

                offset += FL.dyn_bit_pack(delta_buf, out[offset..], width);
            }

            return byte_offset + N_BYTES * offset;
        }

        /// Returns the number of bytes decompressed from the `input`
        fn delta_decompress(
            noalias scratch: *[1024]T,
            noalias transposed: *[1024]T,
            noalias input: []const u8,
            noalias output: []T,
        ) Error!usize {
            if (output.len > std.math.maxInt(u32)) {
                return Error.InvalidInput;
            }
            const len: u32 = @intCast(output.len);

            const n_whole_blocks = len / 1024;
            const n_remainder = len % 1024;

            if (input.len < 1 + n_whole_blocks) {
                return Error.InvalidInput;
            }

            const remainder_width = input[0];
            if (remainder_width > N_BITS) {
                return Error.InvalidInput;
            }

            const block_widths = input[1 .. 1 + n_whole_blocks];

            const remainder_packed_len = (@as(u64, n_remainder) * @as(u64, remainder_width) + N_BITS - 1) / N_BITS;
            var total_packed_len = remainder_packed_len + 1;
            for (block_widths) |w| {
                if (w > N_BITS) {
                    return Error.InvalidInput;
                }
            }
            for (block_widths) |w| {
                total_packed_len += @as(u64, w) * 1024 / N_BITS + FL.N_LANES;
            }

            const data_section_byte_offset = 1 + n_whole_blocks;

            const data_section_bytes = input[data_section_byte_offset..];
            if (data_section_bytes.len < total_packed_len * N_BYTES) {
                return Error.InvalidInput;
            }

            const data_section: []align(1) const U = @ptrCast(data_section_bytes[0 .. total_packed_len * N_BYTES]);

            // number of T read so far, NOT number of bytes
            var offset: usize = 0;

            // read remainder data
            if (n_remainder > 0) {
                const base = data_section[0];
                offset += 1;

                if (IS_SIGNED) {
                    const zigzagged: []U = @ptrCast(scratch[0..n_remainder]);
                    const n_read = try SPack.delta_unpack(
                        base,
                        data_section[offset .. offset + remainder_packed_len],
                        zigzagged,
                        remainder_width,
                    );
                    std.debug.assert(n_read == remainder_packed_len);
                    ZigZag(T).decode(zigzagged, output[0..n_remainder]);
                } else {
                    const scratch_u: []U = @ptrCast(scratch[0..n_remainder]);
                    const n_read = try SPack.delta_unpack(
                        base,
                        data_section[offset .. offset + remainder_packed_len],
                        scratch_u,
                        remainder_width,
                    );
                    std.debug.assert(n_read == remainder_packed_len);
                    @memcpy(output[0..n_remainder], scratch_u);
                }
                offset += remainder_packed_len;
            } else {
                offset += 1;
            }

            // Read whole blocks
            for (0..n_whole_blocks) |block_idx| {
                const width = block_widths[block_idx];

                const bases_p = data_section[offset..];
                const bases: *align(1) const [FL.N_LANES]U = bases_p[0..FL.N_LANES];
                offset += FL.N_LANES;

                const transposed_buf: *[1024]U = @ptrCast(transposed);

                const out_offset = output[block_idx * 1024 + n_remainder ..];
                if (IS_SIGNED) {
                    const zigzagged: *[1024]U = @ptrCast(scratch);
                    offset += FL.dyn_undelta_pack(data_section[offset..], bases, transposed_buf, width);
                    FL.untranspose(transposed_buf, zigzagged);
                    ZigZag(T).decode1024(zigzagged, out_offset[0..1024]);
                } else {
                    const scratch_u: *[1024]U = @ptrCast(scratch);
                    offset += FL.dyn_undelta_pack(data_section[offset..], bases, transposed_buf, width);
                    FL.untranspose(transposed_buf, scratch_u);
                    out_offset[0..1024].* = scratch_u.*;
                }
            }

            std.debug.assert(offset == total_packed_len);
            std.debug.assert(offset == data_section.len);

            return data_section_byte_offset + @sizeOf(T) * offset;
        }

        fn max1024(input: *const [1024]U) U {
            var m = input[0];
            for (0..1024) |i| {
                m = @max(input[i], m);
            }
            return m;
        }

        fn minmax1024(input: *const [1024]U) struct { U, U } {
            var min = input[0];
            var max = input[0];

            for (0..1024) |i| {
                min = @min(input[i], min);
                max = @max(input[i], max);
            }

            return .{ min, max };
        }

        fn needed_width(range: U) u7 {
            return N_BITS - @clz(range);
        }

        /// Load input data, apply zigzag encoding if needed
        /// returns the loaded data and a bit set to 1 if zigzag encoding is applied
        fn load_remainder(input: []const T, scratch: []T) []const U {
            std.debug.assert(input.len == scratch.len);
            if (IS_SIGNED) {
                const zigzagged: []U = @ptrCast(scratch);
                ZigZag(T).encode(input, zigzagged);
                return zigzagged;
            } else {
                return @ptrCast(input);
            }
        }

        /// Load input data, apply zigzag encoding if needed
        /// returns the loaded data and a bit set to 1 if zigzag encoding is applied
        fn load_block(input: *const [1024]T, scratch: *[1024]T) *const [1024]U {
            if (IS_SIGNED) {
                const zigzagged: *[1024]U = @ptrCast(scratch);
                ZigZag(T).encode(input, zigzagged);
                return zigzagged;
            } else {
                return @ptrCast(input);
            }
        }
    };
}

fn Test(comptime T: type) type {
    return struct {
        const MAX_NUM_INTS = 123321;

        const Z = Zint(T);

        const TestCtx = struct {
            const page_allocator = std.heap.page_allocator;

            ctx: Ctx,
            input: []T,
            compressed: []u8,
            output: []T,

            fn compress_bound() usize {
                return @max(
                    Z.bitpack_compress_bound(MAX_NUM_INTS),
                    Z.forpack_compress_bound(MAX_NUM_INTS),
                    Z.deltapack_compress_bound(MAX_NUM_INTS),
                );
            }

            fn init() TestCtx {
                const input = page_allocator.alloc(T, MAX_NUM_INTS) catch unreachable;
                const output = page_allocator.alloc(T, MAX_NUM_INTS) catch unreachable;
                const compressed = page_allocator.alloc(
                    u8,
                    compress_bound(),
                ) catch unreachable;
                const ctx = Ctx.init(page_allocator) catch unreachable;

                return .{
                    .input = input,
                    .output = output,
                    .compressed = compressed,
                    .ctx = ctx,
                };
            }

            fn deinit(self: *TestCtx) void {
                page_allocator.free(self.input);
                page_allocator.free(self.output);
                page_allocator.free(self.compressed);

                self.input = &.{};
                self.output = &.{};
                self.compressed = &.{};

                self.ctx.deinit(page_allocator);
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

        fn roundtrip_bitpack(ctx: TestCtx, input: []const u8) anyerror!void {
            const in = read_input(ctx.input, input);

            const compressed_len = try Z.bitpack_compress(
                ctx.ctx,
                in,
                ctx.compressed,
            );

            std.debug.assert(
                compressed_len <= Z.bitpack_compress_bound(MAX_NUM_INTS),
            );

            const n_decompress = try Z.bitpack_decompress(
                ctx.ctx,
                ctx.compressed[0..compressed_len],
                ctx.output[0..in.len],
            );

            try std.testing.expectEqual(n_decompress, compressed_len);

            try std.testing.expectEqualSlices(T, in, ctx.output[0..in.len]);
        }

        test roundtrip_bitpack {
            var ctx = TestCtx.init();
            defer ctx.deinit();
            try std.testing.fuzz(ctx, roundtrip_bitpack, .{});
        }

        fn roundtrip_forpack(ctx: TestCtx, input: []const u8) anyerror!void {
            const in = read_input(ctx.input, input);

            const compressed_len = try Z.forpack_compress(
                ctx.ctx,
                in,
                ctx.compressed,
            );

            std.debug.assert(
                compressed_len <= Z.forpack_compress_bound(MAX_NUM_INTS),
            );

            const n_decompress = try Z.forpack_decompress(
                ctx.ctx,
                ctx.compressed[0..compressed_len],
                ctx.output[0..in.len],
            );

            try std.testing.expectEqual(n_decompress, compressed_len);

            try std.testing.expectEqualSlices(T, in, ctx.output[0..in.len]);
        }

        test roundtrip_forpack {
            var ctx = TestCtx.init();
            defer ctx.deinit();
            try std.testing.fuzz(ctx, roundtrip_forpack, .{});
        }

        fn roundtrip_delta_pack(ctx: TestCtx, input: []const u8) anyerror!void {
            const in = read_input(ctx.input, input);

            const compressed_len = try Z.deltapack_compress(
                ctx.ctx,
                in,
                ctx.compressed,
            );

            std.debug.assert(
                compressed_len <= Z.deltapack_compress_bound(MAX_NUM_INTS),
            );

            const n_decompress = try Z.deltapack_decompress(
                ctx.ctx,
                ctx.compressed[0..compressed_len],
                ctx.output[0..in.len],
            );

            try std.testing.expectEqual(n_decompress, compressed_len);

            try std.testing.expectEqualSlices(T, in, ctx.output[0..in.len]);
        }

        test roundtrip_delta_pack {
            var ctx = TestCtx.init();
            defer ctx.deinit();
            try std.testing.fuzz(ctx, roundtrip_delta_pack, .{});
        }

        fn garbage_bitpack(ctx: TestCtx, input: []const u8) anyerror!void {
            if (input.len < 4) {
                _ = Z.bitpack_decompress(ctx.ctx, input, ctx.output) catch return;
            } else {
                const n_ints: u32 = @bitCast(@as([*]const [4]u8, @ptrCast(input.ptr))[0]);
                const num_ints = @min(@as(usize, n_ints), MAX_NUM_INTS);
                _ = Z.bitpack_decompress(ctx.ctx, input[4..], ctx.output[0..num_ints]) catch return;
            }
        }

        test garbage_bitpack {
            var ctx = TestCtx.init();
            defer ctx.deinit();
            try std.testing.fuzz(ctx, garbage_bitpack, .{});
        }

        fn garbage_forpack(ctx: TestCtx, input: []const u8) anyerror!void {
            if (input.len < 4) {
                _ = Z.forpack_decompress(ctx.ctx, input, ctx.output) catch return;
            } else {
                const n_ints: u32 = @bitCast(@as([*]const [4]u8, @ptrCast(input.ptr))[0]);
                const num_ints = @min(@as(usize, n_ints), MAX_NUM_INTS);
                _ = Z.forpack_decompress(ctx.ctx, input[4..], ctx.output[0..num_ints]) catch return;
            }
        }

        test garbage_forpack {
            var ctx = TestCtx.init();
            defer ctx.deinit();
            try std.testing.fuzz(ctx, garbage_forpack, .{});
        }

        fn garbage_delta_pack(ctx: TestCtx, input: []const u8) anyerror!void {
            if (input.len < 4) {
                _ = Z.deltapack_decompress(ctx.ctx, input, ctx.output) catch return;
            } else {
                const n_ints: u32 = @bitCast(@as([*]const [4]u8, @ptrCast(input.ptr))[0]);
                const num_ints = @min(@as(usize, n_ints), MAX_NUM_INTS);
                _ = Z.deltapack_decompress(ctx.ctx, input[4..], ctx.output[0..num_ints]) catch return;
            }
        }

        test garbage_delta_pack {
            var ctx = TestCtx.init();
            defer ctx.deinit();
            try std.testing.fuzz(ctx, garbage_delta_pack, .{});
        }
    };
}

test {
    _ = Test(u8);
    _ = Test(u16);
    _ = Test(u32);
    _ = Test(u64);

    _ = Test(i8);
    _ = Test(i16);
    _ = Test(i32);
    _ = Test(i64);

    _ = Test(i128);
    _ = Test(u128);

    _ = Test(i256);
    _ = Test(u256);
}
