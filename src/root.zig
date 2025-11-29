const fastlanes = @import("./fastlanes.zig");
const zint = @import("./zint.zig");

pub const Zint = zint.Zint;

test {
    _ = fastlanes;
    _ = zint;
}
