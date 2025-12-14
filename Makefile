fuzz:
	zig build fuzz --fuzz -Doptimize=ReleaseSafe -j128
fuzz_debug:
	zig build fuzz --fuzz -Doptimize=Debug -j128
bench:
	zig build bench -Doptimize=ReleaseFast -Dcpu=native-avx512f
bench_safe:
	zig build bench -Doptimize=ReleaseSafe -Dcpu=native-avx512f
