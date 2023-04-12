paramer : src/main.cpp src/cmd.cpp src/utils.cpp src/fastx.cpp src/seq.cpp src/kraken2.cpp src/bloom.cpp
	g++ -std=c++17 \
		-I third_party -I third_party/ntHash-2.2.0 -I src -I third_party/cxxopts/include \
		$^ \
		-o $@ -lz -O3

test : tests/unit_tests/seq_tests.cpp tests/unit_tests/fastx_tests.cpp tests/unit_tests/utils_tests.cpp tests/unit_tests/bloom_tests.cpp
	g++ -std=c++17 \
		-I third_party -I third_party/doctest/doctest -I third_party/ntHash-2.2.0 -I src -I third_party/cxxopts/include \
		$^ src/utils.cpp src/fastx.cpp src/seq.cpp src/kraken2.cpp src/bloom.cpp \
		-o $@ -lz
	./$@; rm $@

stress-test : tests/stress_tests/bloom_tests.cpp
	g++ -std=c++17 \
		-I third_party -I third_party/doctest/doctest -I third_party/ntHash-2.2.0 -I src -I third_party/cxxopts/include \
		$^ src/utils.cpp src/fastx.cpp src/seq.cpp src/kraken2.cpp src/bloom.cpp \
		-o $@ -lz -O3
	./$@; rm $@
