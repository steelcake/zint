const fastlanes = @import("./fastlanes.zig");
pub const FastLanes = fastlanes.FastLanes;

const zint = @import("./zint.zig");
pub const Zint = zint.Zint;

test {
    _ = fastlanes;
    _ = zint;
}
