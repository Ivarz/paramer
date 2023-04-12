#include <iostream>
#include "cmd.h"


void print_cmd_usage() {
  std::cerr << "Usage:\n";
  std::cerr << "  paramer <COMMAND>\n";
  std::cerr << '\n';
  std::cerr << "  where command is:\n";
  std::cerr << "      mask\t\tmask fasta file with kraken2 file and/or kmers found in another fasta file\n";
  std::cerr << "      bloom-build\tbuild Bloom's filter\n";
  std::cerr << "      bloom-search\tsearch in Bloom's filter\n";
  std::cerr << "      stats\t\tPerform various statistics operations\n";
}

int main(int argc, char **argv) {

  if (argc < 2) {
    print_cmd_usage();
    return 0;
  }

  if (std::string(argv[1]) == "mask") { Cmd::Mask::run(argc, argv); }
  else if (std::string(argv[1]) == "bloom-build") { Cmd::BloomBuild::run(argc, argv); }
  else if (std::string(argv[1]) == "bloom-search") { Cmd::BloomSearch::run(argc, argv); }
  else if (std::string(argv[1]) == "extend") { Cmd::Extend::run(argc, argv); }
  else if (std::string(argv[1]) == "stats") { Cmd::Stats::run(argc, argv); }
  else { print_cmd_usage(); }

  return 0;
}
