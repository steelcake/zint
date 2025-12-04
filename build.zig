const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const fastlanes = b.addModule("fastlanes", .{
        .root_source_file = b.path("fastlanes/root.zig"),
        .target = target,
        .optimize = optimize,
    });

    const zint = b.addModule("zint", .{
        .root_source_file = b.path("zint/root.zig"),
        .target = target,
        .optimize = optimize,
    });
    zint.addImport("fastlanes", fastlanes);

    const tests = b.createModule(.{
        .root_source_file = b.path("tests.zig"),
        .target = target,
        .optimize = optimize,
    });

    const tests_target = b.addTest(.{
        .root_module = tests,
    });
    const run_tests = b.addRunArtifact(tests_target);

    const test_step = b.step("test", "Run tests");
    test_step.dependOn(&run_tests.step);

    const fuzz_target = b.addTest(.{
        .root_module = tests,
        // Required for running fuzz tests
        // https://github.com/ziglang/zig/issues/23423
        .use_llvm = true,
    });
    const run_fuzz = b.addRunArtifact(fuzz_target);

    const fuzz_step = b.step("fuzz", "Run fuzz tests");
    fuzz_step.dependOn(&run_fuzz.step);
}
