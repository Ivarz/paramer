#include "cmd.h"
#include "bloom.h"
#include "fastx.h"
#include "kraken2.h"
#include "seq.h"
#include "utils.h"
#include <cxxopts.hpp>
#include <iostream>
#include <vector>

namespace Cmd {
void print_help(const cxxopts::Options &options) {
  std::cerr << options.help();
}
namespace Mask {
	int run(int argc, char **argv) {

	  cxxopts::Options options(
		  "mask",
		  "Mask fasta file with a corresponding kraken2 classification file");
	  options.add_options()("f,fasta", "Fasta file to mask",
							cxxopts::value<std::string>())(
		  "k,kraken2",
		  "Kraken2 to use for masking. Can be supplied multiple times.",
		  cxxopts::value<std::vector<std::string>>())(
		  "r,reference",
		  "Reference file to use for masking in fasta format. Can be supplied "
		  "multiple times.",
		  cxxopts::value<std::vector<std::string>>())(
		  "klen", "Kmer length to use for filter",
		  cxxopts::value<size_t>()->default_value("31"))("h,help",
														   "Help message");

	  if (argc < 3) {
		print_help(options);
		return 0;
	  }

	  auto result = options.parse(argc - 1, argv + 1);
	  if (result.count("help")) {
		std::cerr << "Showing help message\n";
		print_help(options);
		return 0;
	  }

	  std::string fasta_fname = result["fasta"].as<std::string>();
	  std::vector<std::string> kraken2_fnames =
		  result["kraken2"].as<std::vector<std::string>>();
	  std::vector<std::string> reference_fnames =
		  result["reference"].as<std::vector<std::string>>();

	  size_t kmer_size = result["klen"].as<size_t>();

	  if (fasta_fname.size() == 0) {
		std::cerr << "Provide fasta file for masking\n";
		print_help(options);
		return 1;
	  }

	  Fasta::loadSoftmaskAndPrint(fasta_fname, kraken2_fnames, reference_fnames, kmer_size);
	  return 0;
	}

} // namespace Mask

namespace BloomBuild {
	int run(int argc, char **argv) {

	  cxxopts::Options options("bloom-build", "build Bloom's filter");
	  options.add_options()(
		  "r,reference",
		  "Reference in fasta or fastq format. Can be supplied multiple times",
		  cxxopts::value<std::vector<std::string>>())(
		  "k,klen", "Kmer length to use for filter", cxxopts::value<uint64_t>())(
		  "s,size", "filter size in bytes",
		  cxxopts::value<uint64_t>()->default_value("100000000"))(
		  "l,seqlen", "Minimum sequence length to use for inserting in filter",
		  cxxopts::value<uint64_t>()->default_value("100"))(
		  "n,nhash", "Number of hashes to use", cxxopts::value<uint64_t>())(
		  "o,output", "output file", cxxopts::value<std::string>())("h,help",
																	"Help message");

	  if (argc < 3) {
		print_help(options);
		return 0;
	  }

	  auto result = options.parse(argc - 1, argv + 1);
	  if (result.count("help")) {
		std::cerr << "Showing help message\n";
		print_help(options);
		return 0;
	  }

	  std::vector<std::string> seq_fnames =
		  result["reference"].as<std::vector<std::string>>();
	  std::string output_fname = result["output"].as<std::string>();

	  if (seq_fnames.size() == 0) {
		std::cerr << "Provide fasta/fastq input\n";
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
	  // TODO add canonical kmers (lexicographically smaller min(kmer,
	  // revcom(kmer)))
	  for (const std::string &fname : seq_fnames) {
		std::cerr << fname << '\n';
		std::optional<FileFormat> fformat = Fastx::inferFileFormat(fname);
		if (fformat) {
		  switch (*fformat) {
		  case FileFormat::Fasta:
			std::cerr << fname << " inferred as fasta\n";
			blmf.addFasta(fname, seqlen);
			break;
		  case FileFormat::Fastq:
			std::cerr << fname << " inferred as fastq\n";
			blmf.addFastq(fname, seqlen);
			break;
		  default:
			std::cerr << "Unhandled file format\n";
		  }
		} else {
		  std::cerr << "Unrecognized file format\n";
		}
	  }
	  blmf.write(output_fname);

	  return 0;
	}
} // namespace BloomBuild

namespace BloomSearch {
	int run(int argc, char **argv) {

	  cxxopts::Options options("bloom-search", "Search in Bloom's filter");
	  options.add_options()("b,bloom", "Bloom filter",
							cxxopts::value<std::string>())(
		  "s,sequence", "Sequence to search",
		  cxxopts::value<std::string>()->default_value(""))(
		  "1,mates1", "Reads mates 1 (fastq(.gz))", cxxopts::value<std::string>())(
		  "2,mates2", "Reads mates 2 (fastq(.gz))", cxxopts::value<std::string>())(
		  "c,mincount", "minimum number of matching kmers for hit",
		  cxxopts::value<size_t>()->default_value("50"))(
		  "u,unpaired", "Single end reads",
		  cxxopts::value<std::string>())("h,help", "Help message");

	  if (argc < 3) {
		print_help(options);
		return 0;
	  }
	  auto result = options.parse(argc - 1, argv + 1);
	  std::string seq = result["sequence"].as<std::string>();
	  std::string bloom_filter_name = result["bloom"].as<std::string>();

	  std::optional<Bloom::Filter> bloom_filter = Bloom::Filter::load(bloom_filter_name);
	  size_t hits = 0;
	  if (seq.size() > 0) {
		hits = bloom_filter->searchSeq(seq);
		std::cout << "Kmers matching: " << hits << '\n';
		return 0;
	  }

	  std::string mates1_fname = result["mates1"].as<std::string>();
	  std::string mates2_fname = result["mates2"].as<std::string>();
	  size_t hit_threshold = result["mincount"].as<size_t>();

	  Gz::Reader mates1_reader = Gz::Reader(mates1_fname);
	  Gz::Reader mates2_reader = Gz::Reader(mates2_fname);

	  std::optional<Fastq::Pair> curr_rec_pair =
		  Fastq::nextRecordPair(mates1_reader, mates2_reader);

	  while (curr_rec_pair) {
		size_t kmer_hits = bloom_filter->searchFastqPair(*curr_rec_pair);
		if (kmer_hits >= hit_threshold) {
		  // std::cout << "Kmer hits\t" << kmer_hits << '\n';
		  curr_rec_pair->first.print();
		  curr_rec_pair->second.print();
		}
		curr_rec_pair = Fastq::nextRecordPair(mates1_reader, mates2_reader);
	  }

	  return 0;
	}
} // namespace BloomSearch

namespace Extend {
	int run(int argc, char **argv) {

	  cxxopts::Options options("extend",
							   "Extend sequence in 3' and 5' directions with kmers found in Bloom's filter");
	  options.add_options()("b,bloom", "Bloom filter",
							cxxopts::value<std::string>())(
		  "s,sequence", "Sequence to search",
		  cxxopts::value<std::string>()->default_value(""))("h,help",
															"Help message");

	  if (argc < 3) {
		print_help(options);
		return 0;
	  }
	  auto result = options.parse(argc - 1, argv + 1);
	  std::string seq = result["sequence"].as<std::string>();
	  std::string bloom_filter_name = result["bloom"].as<std::string>();

	  std::optional<Bloom::Filter> bloom_filter = Bloom::Filter::load(bloom_filter_name);
	  if (bloom_filter) {
		  std::vector<std::string> candidate_seqs =  bloom_filter->extendSeq(seq);
		  for (const auto& seq: candidate_seqs) {
			  std::cout << seq << '\n';
		  }
	  } else {
		  std::cerr << "Failed to load bloom filter\n";
	  }

	  return 0;
	}
} // namespace Extend

namespace Stats {
	int run(int argc, char **argv) {

		cxxopts::Options options("stats",
							   "Perform various statistics operations");
			options.add_options()("b,bloom", "Bloom filter", cxxopts::value<std::string>())
							   ("h,help", "Help message");

			if (argc < 3) {
			print_help(options);
				return 0;
			}
			auto result = options.parse(argc - 1, argv + 1);
			std::string bloom_filter_name = result["bloom"].as<std::string>();
			std::optional<Bloom::Filter> bloom_filter = Bloom::Filter::load(bloom_filter_name);
			std::cout << "bloom filter size:\t" << bloom_filter->size() << '\n';
			std::cout << "set bits:\t" << bloom_filter->setBitsCount() << '\n';
			std::cout << "false positive rate:\t" << bloom_filter->falsePostiveRate() << '\n';

			return 0;
		}
	} // namespace Stats
} // namespace Cmd

