paramer : src/main.cpp src/utils.cpp src/fastx.cpp src/seq.cpp src/kraken2.cpp
	g++ -std=c++17 -I third_party -I src -I third_party/cxxopts/include $< src/utils.cpp src/fastx.cpp src/seq.cpp src/kraken2.cpp -o $@ -lz
