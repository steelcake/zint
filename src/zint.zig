const ZigZag = @import("./zigzag.zig").ZigZag;
const FastLanes = @import("./fastlanes.zig").FastLanes;

pub fn Zint(comptime T: type) type {
    _ = T;
    return struct {};
}

fn Signed(comptime T: type) type {
    switch (T) {
        i8, i16, i32, i64 => {},
        else => @compileError("unexpected type"),
    }

    // const ZZ = ZigZag(T);

    return struct {

    };
}

fn Impl(comptime T: type) type {
    switch (T) {
        u8, u16, u32, u64 => {},
        else => @compileError("unexpected type."),
    }

    const FL = FastLanes(T);

    return struct {
        pub fn bitpack_compress(noalias input: []const T, noalias output: []u8) void {
            var offset: usize = 0;

            const n_whole_blocks = input.len / 1024;
            const whole_blocks: []const [1024]T = @ptrCast(input[0 .. n_whole_blocks * 1024]);
            for (whole_blocks) |block| {
                offset += FL.dyn_bit_pack(block, out[offset..]);
            }
        }

        pub fn bitpack_decompress() void {}
        pub fn bitpack_compress_bound() void {}

        pub fn forpack_compress() void {}
        pub fn forpack_decompress() void {}
        pub fn forpack_compress_bound() void {}

        pub fn delta_compress() void {}
        pub fn delta_decompress() void {}
        pub fn delta_compress_bound() void {}
    };
}

pub const Context = struct {
    buf: *[1024]u256,
};
