#include <iostream>
#include <cxxopts.hpp>
#include <vector>
#include "utils.h"
#include "fastx.h"
#include "seq.h"
#include "kraken2.h"

namespace Cmd {
	namespace Mask {

		void print_help(const cxxopts::Options& options) {
			std::cerr << options.help();
		}

		std::vector<Fasta::Rec> softmaskFastaWithKraken2(const std::string& fa_fname, const std::vector<std::string>& k2_fnames) {

			Gz::Reader gzrfa(fa_fname);
			std::vector<Gz::Reader> k2_readers;

			for (const auto& fn: k2_fnames) {
				k2_readers.emplace_back(fn);
			}

			std::optional<Fasta::Rec> fa_rec = Fasta::nextRecord(gzrfa);
			std::vector<std::optional<Kraken2::Rec>> k2_recs;

			Kraken2::nextRecord(k2_readers[0]);

			for (auto& gzkr2: k2_readers) {
				k2_recs.push_back(Kraken2::nextRecord(gzkr2));
			}

			std::vector<Fasta::Rec> results;

			while (fa_rec) {
				for (const auto& k2_rec: k2_recs) {
					if (k2_rec) {
						size_t min_size = fa_rec->seq_id.size() > k2_rec->seq_id.size()
							? k2_rec->seq_id.size()
							: fa_rec->seq_id.size();
						if (fa_rec->seq_id.substr(0, min_size) == k2_rec->seq_id.substr(0, min_size)) {
							fa_rec->softmaskWithKraken2(*k2_rec);
							results.push_back(*fa_rec);
						} else {
							std::cerr << "Ids do not match\n";
							std::cerr << fa_rec->seq_id.substr(min_size) << "\t" << k2_rec->seq_id.substr(min_size) << '\n';
						}
					}
				}

				fa_rec = Fasta::nextRecord(gzrfa);
				k2_recs.clear();

				for (auto& gzkr2: k2_readers) {
					k2_recs.push_back(Kraken2::nextRecord(gzkr2));
				}

			}
			return results;
		}

		int run(int argc, char** argv) {

			cxxopts::Options options("mask", "Mask fasta file with a corresponding kraken2 classification file");
			options.add_options()
				("f,fasta", "Fasta file to mask", cxxopts::value<std::string>())
				("k,kraken2", "Kraken2 to use for masking. Can be supplied multiple times.", cxxopts::value<std::vector<std::string>>())
				("h,help", "Help message");

			if (argc < 3) {
				print_help(options);
				return 0;
			}

			auto result = options.parse(argc-1, argv+1);
			if (result.count("help")) {
				std::cerr << "Showing help message\n";
				print_help(options);
				return 0;
			}

			std::string fasta_fname = result["fasta"].as<std::string>();
			std::vector<std::string> kraken2_fnames = result["kraken2"].as<std::vector<std::string>>();

			if (fasta_fname.size() == 0) {
				std::cerr << "Provide fasta input\n"; 
				print_help(options);
				return 1;
			}
			if (kraken2_fnames.size() == 0) {
				std::cerr << "Provide kraken2 input\n"; 
				print_help(options);
				return 1;
			}

			auto masked_fasta = softmaskFastaWithKraken2(fasta_fname, kraken2_fnames);
			for (const auto& fa_rec: masked_fasta) {
				fa_rec.print();
			}
			return 0;
		}

	}
}

void print_cmd_usage() {
	std::cerr << "Usage:\n";
	std::cerr << "  paramer <COMMAND>\n";
	std::cerr << '\n';
	std::cerr << "  where command is:\n";
	std::cerr << "      mask\tmask fasta file with kraken2 file\n";
}

int main(int argc, char** argv) {

	std::string templ_str = "ATCGTctagtagctgcatgCTCAGCGTCNGACGACGTAcgatcgagCAGTCnATCGATgatgatgagaCGAaA";
	Fasta::Rec fr = Fasta::Rec("test_id", templ_str);
	auto res = fr.splitOnMask();
	for (const auto& r: res) {
		r.print();
	}
	//while (reg.first != templ_str.size()) {
		////std::cout << reg.first << ".." << reg.second << '\n';
		//std::cout << templ_str.substr(reg.first, reg.second - reg.first) << '\n';
		//reg = Dna::nextToggleMaskedRegion(templ_str, reg.second);
	//}
	//for (const auto& s: splitted) {
		//std::cout << s << '\n';
	//}
	return EXIT_SUCCESS;
	if (argc < 2) {
		print_cmd_usage();
		return 0;
	}
	if (std::string(argv[1]) == "mask") {
		Cmd::Mask::run(argc, argv);
		return 0;
	}

	////cxxopts::Options options("paramer", "Foo bar");
	////auto result = options.parse(argc-1, argv+1);
	//std::string file_name = "sample_data/test.unzipped.1.fq";
	//Gz::Reader gzr(file_name);
	//auto fq = Fastq::nextRecord(gzr);

	//while (fq) {
		//fq->print();
		//fq = Fastq::nextRecord(gzr);
	//}
	

	return 0;
}

