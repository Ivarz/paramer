#include <iostream>
#include <cxxopts.hpp>
#include <vector>
#include "utils.h"
#include "fastx.h"
#include "seq.h"
#include "kraken2.h"
#include "bloom.h"

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
	namespace BloomBuild {
		void print_help(const cxxopts::Options& options) {
			std::cerr << options.help();
		}
		int run(int argc, char** argv) {

			cxxopts::Options options("bloom-build", "build Bloom's filter");
			options.add_options()
				("r,reference", "Reference in fasta format", cxxopts::value<std::string>())
				("k,klen", "Kmer length to use for filter", cxxopts::value<uint64_t>())
				("s,size", "filter size in bytes", cxxopts::value<uint64_t>()->default_value("100000000"))
				("l,seqlen", "Minimum sequence length to use for inserting in filter", cxxopts::value<uint64_t>()->default_value("100"))
				("n,nhash", "Number of hashes to use", cxxopts::value<uint64_t>())
				("o,output", "output file", cxxopts::value<std::string>())
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

			std::string fasta_fname = result["reference"].as<std::string>();
			std::string output_fname = result["output"].as<std::string>();

			if (fasta_fname.size() == 0) {
				std::cerr << "Provide fasta input\n"; 
				print_help(options);
				return 1;
			}
			if (output_fname.size() == 0) {
				std::cerr << "Provide output filename\n"; 
				print_help(options);
				return 1;
			}
			uint64_t klen = result["klen"].as<uint64_t>();
			uint64_t size = result["size"].as<uint64_t>();
			uint64_t seqlen = result["seqlen"].as<uint64_t>();
			uint64_t nhash = result["nhash"].as<uint64_t>();
			std::string output = result["output"].as<std::string>();
			Bloom::Filter blmf = Bloom::Filter(size, klen, nhash);
			blmf.addFasta(fasta_fname, seqlen);
			blmf.write(output_fname);
			Bloom::Filter blmf2 = Bloom::Filter(output_fname);

			return 0;
		}
	}

	namespace BloomSearch {
		void print_help(const cxxopts::Options& options) {
			std::cerr << options.help();
		}
		int run(int argc, char** argv) {

			cxxopts::Options options("bloom-search", "Search in Bloom's filter");
			options.add_options()
				("b,bloom", "Bloom filter", cxxopts::value<std::string>())
				("s,sequence", "Kmer length to use for filter", cxxopts::value<std::string>())
				("h,help", "Help message");

			if (argc < 3) {
				print_help(options);
				return 0;
			}
			auto result = options.parse(argc-1, argv+1);
			std::string seq = result["sequence"].as<std::string>();
			std::string bloom_filter_name = result["bloom"].as<std::string>();

			Bloom::Filter bloom_filter = Bloom::Filter(bloom_filter_name);
			size_t hits = bloom_filter.searchSeq(seq);

			std::cout << "Kmers matching: " << hits << '\n';
			return 0;
		}
	}
}

void print_cmd_usage() {
	std::cerr << "Usage:\n";
	std::cerr << "  paramer <COMMAND>\n";
	std::cerr << '\n';
	std::cerr << "  where command is:\n";
	std::cerr << "      mask\t\tmask fasta file with kraken2 file\n";
	std::cerr << "      bloom-build\tbuild Bloom's filter\n";
	std::cerr << "      bloom-search\tsearch in Bloom's filter\n";
}

int main(int argc, char** argv) {

	//std::string templ_str = "ATCGTctagtagctgcatgCTCAGCGTCNGACGACGTAcgatcgagCAGTCnATCGATgatgatgagaCGAaA";
	//Fasta::Rec fr = Fasta::Rec("test_id", templ_str);
	//auto res = fr.splitOnMask();
	//for (const auto& r: res) {
		//r.print();
	//}
	//while (reg.first != templ_str.size()) {
		////std::cout << reg.first << ".." << reg.second << '\n';
		//std::cout << templ_str.substr(reg.first, reg.second - reg.first) << '\n';
		//reg = Dna::nextToggleMaskedRegion(templ_str, reg.second);
	//}
	//for (const auto& s: splitted) {
		//std::cout << s << '\n';
	//}
	//return EXIT_SUCCESS;
	if (argc < 2) {
		print_cmd_usage();
		return 0;
	}

	if (std::string(argv[1]) == "mask") {
		Cmd::Mask::run(argc, argv);
		return 0;
	} else if (std::string(argv[1]) == "bloom-build") {
		Cmd::BloomBuild::run(argc, argv);
		return 0;
	} else if (std::string(argv[1]) == "bloom-search") {
		Cmd::BloomSearch::run(argc, argv);
		return 0;
	} else {
		print_cmd_usage();
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

