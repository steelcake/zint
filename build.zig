const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const zstd = b.dependency("zstd", .{
        .target = target,
        .optimize = optimize,
    });

    const lz4 = b.dependency("lz4", .{
        .target = target,
        .optimize = optimize,
    });

    const zint = b.addModule("zint", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });

    const tests = b.createModule(.{
        .root_source_file = b.path("src/root.zig"),
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

    const bench = b.addExecutable(.{
        .name = "bench",
        .root_module = b.createModule(.{
            .root_source_file = b.path("bench.zig"),
            .target = target,
            .optimize = optimize,
        }),
    });
    bench.root_module.addImport("zint", zint);
    bench.root_module.linkLibrary(zstd.artifact("zstd"));
    bench.root_module.linkLibrary(lz4.artifact("lz4"));
    const run_bench = b.addRunArtifact(bench);
    const bench_step = b.step("bench", "Run benchmarks");
    bench_step.dependOn(&run_bench.step);

    b.installArtifact(bench);
}
