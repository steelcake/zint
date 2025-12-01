const std = @import("std");
const page_allocator = std.heap.page_allocator;

const Zint = @import("zint").Zint;

const MAX_NUM_INTS = 1 << 20;

fn Context(comptime T: type) type {
    return struct {
        const Self = @This();
        const Z = Zint(T);

        input: []T,
        compressed: []u8,
        output: []T,

        pub fn init() Context(T) {
            const input = page_allocator.alloc(T, MAX_NUM_INTS) catch unreachable;
            const output = page_allocator.alloc(T, MAX_NUM_INTS) catch unreachable;
            const compressed = page_allocator.alloc(
                u8,
                Z.compress_bound(MAX_NUM_INTS),
            ) catch unreachable;

            return .{
                .input = input,
                .output = output,
                .compressed = compressed,
            };
        }

        fn deinit(self: *Self) void {
            page_allocator.free(self.input);
            page_allocator.free(self.output);
            page_allocator.free(self.compressed);

            self.input = &.{};
            self.output = &.{};
            self.compressed = &.{};
        }
    };
}

fn read_input(comptime T: type, typed_input: []T, input: []const u8) []T {
    std.debug.assert(typed_input.len == MAX_NUM_INTS);

    const input_num_ints = input.len / @sizeOf(T);
    const num_ints = @min(input_num_ints, MAX_NUM_INTS);
    const num_bytes = @sizeOf(T) * num_ints;
    @memcpy(@as([]u8, @ptrCast(typed_input))[0..num_bytes], input[0..num_bytes]);

    return typed_input[0..num_ints];
}

fn Roundtrip(comptime T: type) type {
    return struct {
        const size = @sizeOf(T);
        const Z = Zint(T);

        pub fn fuzz_one(ctx: Context(T), input: []const u8) anyerror!void {
            const in = read_input(T, ctx.input, input);

            const compressed_len = Z.compress(
                in,
                ctx.compressed,
            );

            std.debug.assert(
                compressed_len <= Z.compress_bound(MAX_NUM_INTS),
            );

            Z.decompress(
                ctx.compressed[0..compressed_len],
                ctx.output[0..in.len],
            ) catch unreachable;

            std.debug.assert(
                std.mem.eql(T, in, ctx.output[0..in.len]),
            );
        }

        pub fn tst() !void {
            var ctx = Context(T).init();
            defer ctx.deinit();
            try std.testing.fuzz(ctx, fuzz_one, .{});
        }
    };
}

fn Garbage(comptime T: type) type {
    return struct {
        const Z = Zint(T);

        pub fn fuzz_one(ctx: []T, input: []const u8) anyerror!void {
            const output = ctx;

            if (input.len < 4) {
                Z.decompress(input, output) catch return;
            } else {
                const n_ints: u32 = @bitCast(@as([*]const [4]u8, @ptrCast(input.ptr))[0]);
                const num_ints = @min(@as(usize, n_ints), MAX_NUM_INTS);
                Z.decompress(input[4..], output[0..num_ints]) catch return;
            }
        }

        pub fn tst() !void {
            var ctx = Context(T).init();
            defer ctx.deinit();
            try std.testing.fuzz(ctx.input, fuzz_one, .{});
        }
    };
}

test "roundtrip u8" {
    try Roundtrip(u8).tst();
}
test "roundtrip u16" {
    try Roundtrip(u16).tst();
}
test "roundtrip u32" {
    try Roundtrip(u32).tst();
}
test "roundtrip u64" {
    try Roundtrip(u64).tst();
}

test "garbage u8" {
    try Garbage(u8).tst();
}
test "garbage u16" {
    try Garbage(u16).tst();
}
test "garbage u32" {
    try Garbage(u32).tst();
}
test "garbage u64" {
    try Garbage(u64).tst();
}

test "roundtrip i8" {
    try Roundtrip(i8).tst();
}
test "roundtrip i16" {
    try Roundtrip(i16).tst();
}
test "roundtrip i32" {
    try Roundtrip(i32).tst();
}
test "roundtrip i64" {
    try Roundtrip(i64).tst();
}

test "garbage i8" {
    try Garbage(i8).tst();
}
test "garbage i16" {
    try Garbage(i16).tst();
}
test "garbage i32" {
    try Garbage(i32).tst();
}
test "garbage i64" {
    try Garbage(i64).tst();
}
