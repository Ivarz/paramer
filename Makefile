paramer : src/main.cpp src/cmd.cpp src/utils.cpp src/fastx.cpp src/seq.cpp src/kraken2.cpp src/bloom.cpp
	g++ -std=c++17 \
		-I third_party -I third_party/ntHash-2.2.0 -I src -I third_party/cxxopts/include -I third_party/robin-hood-hashing/src/include \
		$^ \
		-o $@ -lz -O3

test : tests/unit_tests/seq_tests.cpp tests/unit_tests/fastx_tests.cpp tests/unit_tests/utils_tests.cpp tests/unit_tests/bloom_tests.cpp
	g++ -std=c++17 \
		-I third_party -I third_party/doctest/doctest -I third_party/ntHash-2.2.0 -I src -I third_party/cxxopts/include  -I third_party/robin-hood-hashing/src/include \
		$^ src/utils.cpp src/fastx.cpp src/seq.cpp src/kraken2.cpp src/bloom.cpp \
		-o $@ -lz
	./$@; rm $@

stress-test : tests/stress_tests/bloom_tests.cpp
	g++ -std=c++17 \
		-I third_party -I third_party/doctest/doctest -I third_party/ntHash-2.2.0 -I src -I third_party/cxxopts/include -I third_party/robin-hood-hashing/src/include \
		$^ src/utils.cpp src/fastx.cpp src/seq.cpp src/kraken2.cpp src/bloom.cpp \
		-o $@ -lz -O3
	./$@; rm $@

paramer_prof : src/main.cpp src/cmd.cpp src/utils.cpp src/fastx.cpp src/seq.cpp src/kraken2.cpp src/bloom.cpp
	mkdir -p profiling
	g++ -pg -std=c++17 \
		-I third_party -I third_party/ntHash-2.2.0 -I src -I third_party/cxxopts/include -I third_party/robin-hood-hashing/src/include \
		$^ \
		-o profiling/$@ -lz -O3

profiles : paramer_prof

	profiling/$< mask \
		-f test_data/enterobius_vermicularis.fa.gz \
		-r test_data_large/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
		-k test_data/evm.stdpf.out.txt \
		-k test_data/evm.uhgg.out.txt > /dev/null 2> /dev/null
		#-r test_data_large/hs_1M.fa.gz \
	mv gmon.out profiling/gmon.mask.out
	gprof profiling/$< profiling/gmon.mask.out > profiling/gprof.mask.out


	profiling/$< bloom-build \
		-r evm.m2.fa.gz \
		-k 31 \
		-n 3 \
		-o tmp.blm
	rm tmp.blm
	mv gmon.out profiling/gmon.build.out
	gprof profiling/$< profiling/gmon.build.out > profiling/gprof.build.out


