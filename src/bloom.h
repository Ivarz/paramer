#ifndef BLOOM_H
#define BLOOM_H
#include "fastx.h"
#include "seq.h"
#include <cstdint>
#include <fstream>
#include <ntHashIterator.hpp>
#include <functional>
#include <vector>
#include <zlib.h>
#include "utils.h"

namespace Bloom {
const size_t BITS_IN_BYTE = 8;
class Filter {
public:
  Filter(uint64_t s, uint64_t k, uint64_t h)
      : filter_size(s), kmer_size(k), hash_n(h) {
    bytevec = std::vector<uint8_t>(s, 0);
  }

  void addSeq(const std::string &seq);
  void addFasta(const std::string &fasta_fname, size_t minsize);
  void addFastq(const std::string &fastq_fname, size_t minsize);
  size_t searchSeq(const std::string &seq) const;
  size_t searchFastqPair(const Fastq::Pair &fq_pair) const;

  std::vector<std::string> extendSeq(const std::string& seq) const;

  void writeRaw(const std::string &out_fname) const;
  int write(const std::string &out_fname) const;
  static std::optional<Filter> loadRaw(const std::string &in_fname);
  static std::optional<Filter> load(const std::string &in_fname);
  size_t size() const { return bytevec.size(); }
  uint64_t hashN() const { return hash_n; }
  uint64_t kmerSize() const { return kmer_size; }
  uint8_t& at(size_t idx) { return bytevec.at(idx); };
  uint8_t atCpy(size_t idx) { return bytevec.at(idx); };

  size_t setBitsCount() const;
  double falsePostiveRate() const;

private:
  uint64_t filter_size; // filter size in bytes
  uint64_t kmer_size;
  uint64_t hash_n;
  std::vector<uint8_t> bytevec;
  void dfs(std::string current_seq,
		  std::unordered_set<uint64_t> seen_kmer_hashes,
		  std::vector<std::string>& candidate_seqs,
		  const std::function<std::string(std::string)>& extract_kmer,
		  const std::function<std::string(std::string, char)>& next_seq
		  ) const;

  void dfs5prime(const std::string& current_seq,
		  std::unordered_set<uint64_t> seen_kmer_hashes,
		  std::vector<std::string>& candidate_seqs
		  ) const;
  void dfs3prime(const std::string& current_seq,
		  std::unordered_set<uint64_t> seen_kmer_hashes,
		  std::vector<std::string>& candidate_seqs
		  ) const;

};
} // namespace Bloom
#endif
