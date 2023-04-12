#ifndef CMD_H
#define CMD_H

#include <cxxopts.hpp>

namespace Cmd {

	void print_help(const cxxopts::Options &options);

	namespace Mask { int run(int argc, char **argv); }
	namespace BloomBuild { int run(int argc, char **argv); }
	namespace BloomSearch { int run(int argc, char **argv); }
	namespace Extend { int run(int argc, char **argv); }
	namespace Stats { int run(int argc, char **argv); }

}

#endif
