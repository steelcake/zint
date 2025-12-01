const native_endian = @import("builtin").target.cpu.arch.endian();

comptime {
    if (native_endian != .little) {
        @compileError("borsh-zig only supports little-endian architectures.");
    }
}

const fastlanes = @import("./fastlanes.zig");
const zint = @import("./zint.zig");

pub const Zint = zint.Zint;

test {
    _ = fastlanes;
    _ = zint;
}
