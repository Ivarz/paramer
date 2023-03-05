#include <iostream>
#include <cxxopts.hpp>
#include "utils.h"
#include "fastx.h"
#include "seq.h"
#include "kraken2.h"

namespace Cmd {
	namespace Mask {

		void print_help() {
			std::cerr << "mask -f FASTA -k KRAKEN2 > out.fa\n";
		}

		std::vector<Fasta::Rec> softmaskFastaWithKraken2(const std::string& fa_fname, const std::string& k2_fname) {

			Gz::Reader gzrfa(fa_fname);
			Gz::Reader gzrk2(k2_fname);
			std::optional<Fasta::Rec> fa_rec = Fasta::nextRecord(gzrfa);
			std::optional<Kraken2::Rec> k2_rec = Kraken2::nextRecord(gzrk2);
			std::vector<Fasta::Rec> results;
			while (fa_rec && k2_rec) {
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
				fa_rec = Fasta::nextRecord(gzrfa);
				k2_rec = Kraken2::nextRecord(gzrk2);
			}
			return results;
		}

		int run(int argc, char** argv) {
			if (argc < 3) {
				print_help();
				return 0;
			}
			cxxopts::Options options("mask", "Mask fasta file with a corresponding kraken2 classification file");
			options.add_options()
				("f,fasta", "Fasta file to mask", cxxopts::value<std::string>()->default_value(""))
				("k,kraken2", "Kraken2 to use for masking", cxxopts::value<std::string>()->default_value(""))
				("h,help", "Help message");

			auto result = options.parse(argc-1, argv+1);
			if (result.count("help")) {
				std::cerr << "Showing help message\n";
				print_help();
				return 0;
			}

			std::string fasta_fname = result["fasta"].as<std::string>();
			std::string kraken2_fname = result["kraken2"].as<std::string>();

			if (fasta_fname.size() == 0) {
				std::cerr << "Provide fasta input\n"; 
				print_help();
				return 1;
			}
			if (kraken2_fname.size() == 0) {
				std::cerr << "Provide kraken2 input\n"; 
				print_help();
				return 1;
			}

			std::cout << fasta_fname << '\t' << kraken2_fname << '\n';
			auto masked_fasta = softmaskFastaWithKraken2(fasta_fname, kraken2_fname);
			for (const auto& fa_rec: masked_fasta) {
				fa_rec.print();
			}


			return 0;
		}

	}
}

void print_cmd_usage() {
	std::cout << "paramer COMMAND\n";
	std::cout << "where command is\n";
	std::cout << "mask\n";
}

int main(int argc, char** argv) {

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

