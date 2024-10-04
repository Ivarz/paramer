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
#include <fstream>

namespace Bloom {
const size_t BITS_IN_BYTE = 8;
enum class Compression {
	RAW,
	GZ,
};

class Filter {
public:
  Filter(uint64_t s, uint64_t k, uint64_t w, uint64_t h)
      : filter_size(s), kmer_size(k), window_size(w), hash_n(h) {
    bytevec = std::vector<uint8_t>(s, 0);
  }

  void addSeq(const std::string& seq);
  void addMinimizers(const std::string& seq);
  void addFasta(const std::string& fasta_fname, size_t minsize);
  void addFastq(const std::string& fastq_fname, size_t minsize);
  size_t searchSeq(const std::string& seq);
  size_t seekSeq(const std::string& seq);
  size_t searchMinimizers(const std::string& seq);
  size_t searchFastqPair(const Fastq::Pair& fq_pair);

  std::vector<std::string> extendSeq(const std::string& seq,
									 int max_candidates,
									 int max_path_length
  );
  std::vector<std::string> extendSeqPair(const std::string& seq1,
										 const std::string& seq2,
										 int max_candidates,
										 int max_path_length
  );

  int writeRaw(const std::string &out_fname) const;
  int writeGz(const std::string &out_fname) const;
  int write(const std::string &out_fname, Compression cmpr) const;
  static std::optional<Filter> loadRaw(const std::string& in_fname);
  static std::optional<Filter> loadGz(const std::string& in_fname);
  static std::optional<Filter> load(const std::string& in_fname);

  static std::optional<Filter> loadPointer(const std::string& in_fname);
  void closePointer();

  static Bloom::Compression inferCompression(const std::string& in_fname);
  size_t size() const { return bytevec.size(); }
  uint64_t hashN() const { return hash_n; }
  uint64_t kmerSize() const { return kmer_size; }
  uint64_t windowSize() const { return window_size; }
  uint8_t at(size_t idx) { return bytevec.at(idx); };
  uint8_t& atRef(size_t idx) { return bytevec.at(idx); };
  uint8_t seekAt(size_t idx);
  uint8_t getByteVecVal(size_t idx) { 
	  if (filter_fptr.is_open()) {
		  return this->seekAt(idx);
	  } else {
		  return this->at(idx);
	  }
  };

  size_t setBitsCount() const;
  double falsePostiveRate() const;


private:
  std::ifstream filter_fptr;
  uint64_t filter_size; // filter size in bytes
  uint64_t kmer_size;
  uint64_t window_size;
  uint64_t hash_n;
  std::vector<uint8_t> bytevec;
  void dfs(std::string current_seq,
		  robin_hood::unordered_set<uint64_t> seen_kmer_hashes,
		  std::vector<std::string>& candidate_seqs,
		  const std::function<std::string(std::string)>& extract_kmer,
		  const std::function<std::string(std::string, char)>& next_seq
		  );

  void dfs5prime(const std::string& current_seq,
		  robin_hood::unordered_set<uint64_t> seen_kmer_hashes,
		  std::vector<std::string>& candidate_seqs
		  );
  void dfs3prime(const std::string& current_seq,
		  robin_hood::unordered_set<uint64_t> seen_kmer_hashes,
		  std::vector<std::string>& candidate_seqs
		  );
  std::vector<std::string> bfs(std::string src,
		  robin_hood::unordered_set<uint64_t> seen_kmer_hashes,
		  const std::function<std::string(std::string)>& extract_kmer,
		  const std::function<std::string(std::string, char)>& next_seq,
		  const std::function<std::string(std::string, std::string)>& add_parent_to_path,
		  const std::function<std::string(std::string, uint64_t)>& clip,
		  int max_candidate_limit,
		  int depth_limit
		  );
  std::vector<std::string> bfs5prime(const std::string& src,
		  robin_hood::unordered_set<uint64_t> seen_kmer_hashes,
		  int max_candidate_limit,
		  int depth_limit
		  );
  std::vector<std::string> bfs3prime(const std::string& src,
		  robin_hood::unordered_set<uint64_t> seen_kmer_hashes,
		  int max_candidate_limit,
		  int depth_limit
		  );

};
} // namespace Bloom
#endif
