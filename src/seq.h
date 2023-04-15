#ifndef SEQ_H
#define SEQ_H

#include <string>
#include <vector>
#include <set>
#include <unordered_set>
#include <algorithm>

namespace Dna {
	using SeqInterval = std::pair<size_t, size_t>;

	std::string revcom(const std::string& seq);

	void softmask(std::string& seq, size_t beg, size_t end);

	void unmask(std::string& seq, size_t beg, size_t end);
	//std::set<std::string> kmers(const std::string& seq, size_t size);
	
	void addKmers(const std::string& seq, size_t size, std::unordered_set<std::string>& kmers);

	void addKmerHashes(const std::string& seq, size_t size, size_t hash_n, std::unordered_set<uint64_t>& kmer_hashes);

	std::vector<uint64_t> getHashes(const std::string& seq, size_t kmer_size, size_t hash_n);

	std::vector<uint64_t> getMinimizerHashes(const std::string& seq, size_t kmer_size, size_t hash_n, size_t window_size);
	std::pair<size_t,size_t> nextToggleMaskedRegion(const std::string& seq, size_t beg);

	std::vector<std::string> splitOnMask(const std::string& seq);

	std::vector<std::pair<SeqInterval, std::string>> splitOnMaskWithInterval(const std::string& seq);

	bool hasMasked(const std::string& seq);

	std::string canonicalKmer(const std::string& kmer);

	void softmaskNotInKmerHashes(std::string& seq
			, const std::unordered_set<uint64_t>& kmer_hashes
			, size_t kmer_size);
}

#endif
