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
    };
}

const Ctx8 = Context(u8);
const Ctx16 = Context(u16);
const Ctx32 = Context(u32);
const Ctx64 = Context(u64);

test "roundtrip 8" {
    var ctx = Ctx8.init();
    defer ctx.deinit();
    try std.testing.fuzz(ctx, Roundtrip(u8).fuzz_one, .{});
}

test "roundtrip 16" {
    var ctx = Ctx16.init();
    defer ctx.deinit();
    try std.testing.fuzz(ctx, Roundtrip(u16).fuzz_one, .{});
}

test "roundtrip 32" {
    var ctx = Ctx32.init();
    defer ctx.deinit();
    try std.testing.fuzz(ctx, Roundtrip(u32).fuzz_one, .{});
}

test "roundtrip 64" {
    var ctx = Ctx64.init();
    defer ctx.deinit();
    try std.testing.fuzz(ctx, Roundtrip(u64).fuzz_one, .{});
}

test "garbage 8" {
    var ctx = Ctx8.init();
    defer ctx.deinit();
    try std.testing.fuzz(ctx.input, Garbage(u8).fuzz_one, .{});
}

test "garbage 16" {
    var ctx = Ctx16.init();
    defer ctx.deinit();
    try std.testing.fuzz(ctx.input, Garbage(u16).fuzz_one, .{});
}

test "garbage 32" {
    var ctx = Ctx32.init();
    defer ctx.deinit();
    try std.testing.fuzz(ctx.input, Garbage(u32).fuzz_one, .{});
}

test "garbage 64" {
    var ctx = Ctx64.init();
    defer ctx.deinit();
    try std.testing.fuzz(ctx.input, Garbage(u64).fuzz_one, .{});
}
