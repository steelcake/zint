pub fn ZigZag(comptime T: type) type {
    return struct {
        pub const ALIGNMENT = 64;
        const N_BITS = @sizeOf(T) * 8;
        pub const U = switch (T) {
            i8 => u8,
            i16 => u16,
            i32 => u32,
            i64 => u64,
            else => @compileError("unsupported type"),
        };

        fn negate(x: U) U {
            return ~x +% 1;
        }

        pub fn encode(
            noalias input: []align(ALIGNMENT) const T,
            noalias output: []align(ALIGNMENT) U,
        ) void {
            for (input, output) |*i, *o| {
                const v = i.*;
                o.* = @bitCast((v >> (N_BITS - 1)) ^ (v << 1));
            }
        }

        pub fn decode(
            noalias input: []align(ALIGNMENT) const U,
            noalias output: []align(ALIGNMENT) T,
        ) void {
            for (input, output) |*i, *o| {
                const v = i.*;
                o.* = @bitCast((v >> 1) ^ negate(v & 1));
            }
        }

        pub fn encode1024(
            noalias input: *align(ALIGNMENT) const [1024]T,
            noalias output: *align(ALIGNMENT) [1024]U,
        ) void {
            for (0..1024) |i| {
                const v = input[i];
                output[i] = @bitCast((v >> (N_BITS - 1)) ^ (v << 1));
            }
        }

        pub fn decode1024(
            noalias input: *align(ALIGNMENT) const [1024]U,
            noalias output: *align(ALIGNMENT) [1024]T,
        ) void {
            for (0..1024) |i| {
                const v = input[i];
                output[i] = @bitCast((v >> 1) ^ negate(v & 1));
            }
        }
    };
}

fn Test(comptime T: type) type {
    return struct {
        const std = @import("std");

        const Z = ZigZag(T);

        const ALIGNMENT = 64;

        fn alloc(comptime Y: type, len: usize) []align(ALIGNMENT) Y {
            return std.heap.page_allocator.alignedAlloc(
                Y,
                std.mem.Alignment.fromByteUnits(ALIGNMENT),
                len,
            ) catch unreachable;
        }

        fn read_input1024(input: []const u8) ?[1024]T {
            if (input.len < @sizeOf(T) * 1024) return null;
            return @bitCast(input[0 .. @sizeOf(T) * 1024].*);
        }

        fn fuzz_zigzag1024(_: void, input: []const u8) anyerror!void {
            const in: [1024]T align(ALIGNMENT) = read_input1024(input) orelse return;

            var output: [1024]Z.U align(ALIGNMENT) = std.mem.zeroes([1024]Z.U);
            Z.encode1024(&in, &output);

            var out: [1024]T align(ALIGNMENT) = std.mem.zeroes([1024]T);
            Z.decode1024(&output, &out);

            try std.testing.expect(std.mem.eql(T, &out, &in));
        }

        test fuzz_zigzag1024 {
            try std.testing.fuzz({}, fuzz_zigzag1024, .{});
        }

        const MAX_INPUT = 12332;
        const Ctx = struct {
            input: []align(ALIGNMENT) T,
            zigzagged: []align(ALIGNMENT) Z.U,
            output: []align(ALIGNMENT) T,
        };

        fn read_input(ctx: Ctx, input: []const u8) []align(ALIGNMENT) const T {
            const in_len = input.len / @sizeOf(T) * @sizeOf(T);
            const len = @min(in_len, MAX_INPUT * @sizeOf(T));

            const o: []align(ALIGNMENT) u8 = @ptrCast(ctx.input);

            @memcpy(o[0..len], input[0..len]);

            return @ptrCast(o[0..len]);
        }

        fn fuzz_zigzag(ctx: Ctx, input: []const u8) anyerror!void {
            const in = read_input(ctx, input);

            Z.encode(in, ctx.zigzagged[0..in.len]);
            Z.decode(ctx.zigzagged[0..in.len], ctx.output[0..in.len]);

            try std.testing.expect(std.mem.eql(T, ctx.output[0..in.len], in));
        }

        test fuzz_zigzag {
            const ctx = Ctx{
                .input = alloc(T, MAX_INPUT),
                .zigzagged = alloc(Z.U, MAX_INPUT),
                .output = alloc(T, MAX_INPUT),
            };
            defer std.heap.page_allocator.free(ctx.zigzagged);
            defer std.heap.page_allocator.free(ctx.output);
            try std.testing.fuzz(ctx, fuzz_zigzag, .{});
        }
    };
}

test "i8" {
    _ = Test(i8);
}

test "i16" {
    _ = Test(i16);
}

test "i32" {
    _ = Test(i32);
}

test "i64" {
    _ = Test(i64);
}
