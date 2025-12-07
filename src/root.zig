const native_endian = @import("builtin").target.cpu.arch.endian();

comptime {
    if (native_endian != .little) {
        @compileError("zint only supports little-endian architectures.");
    }
}

const fastlanes = @import("./fastlanes.zig");
const zigzag = @import("./zigzag.zig");
const scalar_bitpack = @import("./scalar_bitpack.zig");
// const zint = @import("./zint.zig");

pub const FastLanes = fastlanes.FastLanes;
// pub const Zint = zint.Zint;

test {
    _ = fastlanes;
    _ = zigzag;
    _ = scalar_bitpack;
    // _ = zint;
}
