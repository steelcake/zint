const fastlanes = @import("./fastlanes.zig");
pub const FastLanes = fastlanes.FastLanes;

const zigzag = @import("./zigzag.zig");
pub const ZigZag = zigzag.ZigZag;

const zint = @import("./zint.zig");
pub const Zint = zint.Zint;

test {
    _ = fastlanes;
    _ = zigzag;
    _ = zint;
}
