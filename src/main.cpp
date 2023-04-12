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
void softmaskFastaWithKraken2(const std::string &fa_fname,
                              const std::vector<std::string> &k2_fnames) {

  Gz::Reader gzrfa(fa_fname);
  std::vector<Gz::Reader> k2_readers;

  for (const auto &fn : k2_fnames) {
    k2_readers.emplace_back(fn);
  }

  std::optional<Fasta::Rec> fa_rec = Fasta::nextRecord(gzrfa);
  std::vector<std::optional<Kraken2::Rec>> k2_recs;

  for (auto &gzkr2 : k2_readers) {
    k2_recs.push_back(Kraken2::nextRecord(gzkr2));
  }

  std::vector<Fasta::Rec> results;

  while (fa_rec) {
    std::cerr << "Processing " << fa_rec->seq_id << '\n';
    for (const auto &k2_rec : k2_recs) {
      if (k2_rec) {
        size_t min_size = fa_rec->seq_id.size() > k2_rec->seq_id.size()
                              ? k2_rec->seq_id.size()
                              : fa_rec->seq_id.size();
        if (fa_rec->seq_id.substr(0, min_size) ==
            k2_rec->seq_id.substr(0, min_size)) {
          std::cerr << "Processing " << fa_rec->seq_id << '\n';
          fa_rec->softmaskWithKraken2(*k2_rec);
          // results.push_back(*fa_rec);
        } else {
          std::cerr << "Ids do not match\n";
          std::cerr << fa_rec->seq_id.substr(0, min_size) << "\t"
                    << k2_rec->seq_id.substr(0, min_size) << '\n';
        }
      }
    }

    fa_rec->print();
    fa_rec = Fasta::nextRecord(gzrfa);
    k2_recs.clear();

    for (auto &gzkr2 : k2_readers) {
      k2_recs.push_back(Kraken2::nextRecord(gzkr2));
    }
  }
  // return results;
}

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
  //softmaskFastaWithKraken2(fasta_fname, kraken2_fnames);
  // auto masked_fasta = softmaskFastaWithKraken2(fasta_fname, kraken2_fnames);
  // for (const auto& fa_rec: masked_fasta) {
  // 	fa_rec.print();
  // }
  //
  //std::unordered_set<uint64_t> fa_kmers =
      //Fasta::loadUnmaskedKmerHashes(fasta_fname, kmer_size);
  //std::cerr << "kmer hashes loaded: " << fa_kmers.size() << '\n';

  //std::cerr << "Dropping hashes found in refs\n";
  //for (std::string ref_name : reference_fnames) {
    //Fasta::dropKmerHashesFound(ref_name, kmer_size, fa_kmers);
  //}
  //std::cerr << "kmer hashes kept: " << fa_kmers.size() << '\n';
  //std::cerr << "clean target fasta file\n";
  //std::vector<Fasta::Rec> masked =
      //Fasta::softmaskNotInKmerHashes(fasta_fname, kmer_size, fa_kmers);
  //for (auto rec : masked) {
    //rec.print();
  //}
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

  cxxopts::Options options("reconstruct",
                           "Reconstruct from seed and Bloom's filter");
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
  size_t hits = 0;
  if (seq.size() > 0) {
    hits = bloom_filter->searchSeq(seq);
    std::cout << "Kmers matching: " << hits << '\n';
  }
  size_t prev_hits = 0;
  std::string new_seq = "";
  std::string prev_seq = "";

  while (seq.size() > prev_seq.size()) {
    prev_seq = seq;
    for (char c : "ATGC") {
      new_seq = seq + c;
      std::cerr << "trying" << '\t' << new_seq << '\n';
      size_t new_hits = bloom_filter->searchSeq(new_seq);
      if (new_hits > hits) {
        seq = new_seq;
        hits = new_hits;
        std::cerr << new_hits << '\t' << seq << '\n';
        continue;
      }
    }
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
