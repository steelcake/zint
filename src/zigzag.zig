pub fn ZigZag(comptime T: type) type {
    return struct {
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
            noalias input: []const T,
            noalias output: []U,
        ) void {
            for (input, output) |*i, *o| {
                const v = i.*;
                o.* = @bitCast((v >> (N_BITS - 1)) ^ (v << 1));
            }
        }

        pub fn decode(
            noalias input: []const U,
            noalias output: []T,
        ) void {
            for (input, output) |*i, *o| {
                const v = i.*;
                o.* = @bitCast((v >> 1) ^ negate(v & 1));
            }
        }

        pub fn encode1024(
            noalias input: *const [1024]T,
            noalias output: *[1024]U,
        ) void {
            for (0..1024) |i| {
                const v = input[i];
                output[i] = @bitCast((v >> (N_BITS - 1)) ^ (v << 1));
            }
        }

        pub fn decode1024(
            noalias input: *const [1024]U,
            noalias output: *[1024]T,
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

        fn read_input1024(input: []const u8) ?[1024]T {
            if (input.len < @sizeOf(T) * 1024) return null;
            return @bitCast(input[0 .. @sizeOf(T) * 1024].*);
        }

        fn fuzz_zigzag1024(_: void, input: []const u8) anyerror!void {
            const in = read_input1024(input) orelse return;

            var output = std.mem.zeroes([1024]Z.U);
            Z.encode1024(&in, &output);

            var out = std.mem.zeroes([1024]T);
            Z.decode1024(&output, &out);

            try std.testing.expect(std.mem.eql(T, &out, &in));
        }

        test fuzz_zigzag1024 {
            try std.testing.fuzz({}, fuzz_zigzag1024, .{});
        }

        const MAX_INPUT = 12332;
        const Ctx = struct {
            input: []T,
            zigzagged: []Z.U,
            output: []T,
        };

        fn read_input(ctx: Ctx, input: []const u8) []const T {
            const in_len = input.len / @sizeOf(T) * @sizeOf(T);
            const len = @min(in_len, MAX_INPUT * @sizeOf(T));

            const o: []align(@sizeOf(T)) u8 = @ptrCast(ctx.input);

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
                .input = try std.heap.page_allocator.alloc(T, MAX_INPUT),
                .zigzagged = try std.heap.page_allocator.alloc(Z.U, MAX_INPUT),
                .output = try std.heap.page_allocator.alloc(T, MAX_INPUT),
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
