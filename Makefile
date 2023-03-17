paramer : src/main.cpp src/utils.cpp src/fastx.cpp src/seq.cpp src/kraken2.cpp src/bloom.cpp
	g++ -std=c++17 \
		-I third_party -I third_party/ntHash-2.2.0 -I src -I third_party/cxxopts/include \
		$< src/utils.cpp src/fastx.cpp src/seq.cpp src/kraken2.cpp src/bloom.cpp \
		-o $@ -lz
