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
            noalias input: *const [1024]T,
            noalias output: *[1024]U,
        ) void {
            for (0..1024) |i| {
                const v = input[i];
                output[i] = @bitCast((v >> (N_BITS - 1)) ^ (v << 1));
            }
        }

        pub fn decode(
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

        fn read_input(input: []const u8) ?[1024]T {
            if (input.len < @sizeOf(T) * 1024) return null;
            return @bitCast(input[0 .. @sizeOf(T) * 1024].*);
        }

        fn fuzz_zigzag(_: void, input: []const u8) anyerror!void {
            const in = read_input(input) orelse return;

            var output = std.mem.zeroes([1024]Z.U);
            Z.encode(&in, &output);

            var out = std.mem.zeroes([1024]T);
            Z.decode(&output, &out);

            try std.testing.expect(std.mem.eql(T, &out, &in));
        }

        test "fuzz_zigzag" {
            try std.testing.fuzz({}, fuzz_zigzag, .{});
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
