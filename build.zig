const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const lib_mod = b.addModule("zint", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });

    const lib_unit_tests = b.addTest(.{
        .root_module = lib_mod,
    });
    const run_lib_unit_tests = b.addRunArtifact(lib_unit_tests);
    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_lib_unit_tests.step);

    const fuzz = b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/fuzz.zig"),
            .target = target,
            .optimize = optimize,
        }),
        // Required for running fuzz tests
        // https://github.com/ziglang/zig/issues/23423
        .use_llvm = true,
    });
    fuzz.root_module.addImport("zint", lib_mod);

    const run_fuzz = b.addRunArtifact(fuzz);

    const fuzz_step = b.step("fuzz", "run fuzz tests");
    fuzz_step.dependOn(&run_fuzz.step);
}
