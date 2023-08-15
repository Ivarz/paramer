#include "seq.h"
#include <algorithm>
#include <iostream>
#include "ntHashIterator.hpp"
#include <limits>
#include <queue>
#include <map>
#include <math.h>

namespace Dna {
	std::string revcom(const std::string& seq) {
		std::string result(seq.size(), 'N');
		if (seq.size() == 0) {
			return result;
		}
		for (size_t i = 0; i < seq.size(); i++) {
			int res_idx = seq.size()-1 - i;
			switch (seq[i]) {
				case 'A':
					result[res_idx] = 'T';
					break;
				case 'a':
					result[res_idx] = 't';
					break;
				case 'T':
					result[res_idx] = 'A';
					break;
				case 't':
					result[res_idx] = 'a';
					break;
				case 'G':
					result[res_idx] = 'C';
					break;
				case 'g':
					result[res_idx] = 'c';
					break;
				case 'C':
					result[res_idx] = 'G';
					break;
				case 'c':
					result[res_idx] = 'g';
					break;
				default:
					result[res_idx] = seq[i];
			}
		}
		return result;
	}

	void softmask(std::string& seq, size_t beg, size_t end) {
		for (size_t i=beg; i<end; i++){
			if (seq[i] < 97) {
				seq[i] += 32;
			}
		}
	}

	void unmask(std::string& seq, size_t beg, size_t end) {
		for (size_t i=beg; i<end; i++){
			if (seq[i] >= 97) {
				seq[i] -= 32;
			}
		}
	}

	void addKmers(const std::string& seq, size_t size, robin_hood::unordered_set<std::string>& kmers) {
		if (seq.size() < size) {
			return; 
		}
		size_t beg = 0;
		size_t end = beg+size;
		while (end < seq.size()) {
			std::string kmer = seq.substr(beg, size);
			kmers.insert(kmer);
			beg++;
			end = beg+size;
		}
	}

	void addKmerHashes(const std::string& seq,
			size_t size,
			size_t hash_n,
			robin_hood::unordered_set<uint64_t>& kmer_hashes
			){
		for (const auto& hash: getHashes(seq, size, hash_n)){
			kmer_hashes.insert(hash);
		}
	}
	bool isMasked(char c) {
		return c == 'a' || c == 't' || c == 'g' || c == 'c' || c == 'N' || c == 'n';

	}
	std::pair<size_t,size_t> nextToggleMaskedRegion(const std::string& seq, size_t beg) {
		if (beg >= seq.size()) {
			return std::pair<size_t,size_t>(seq.size(), seq.size());
		}
		char c = seq[beg];
		bool curr_mask = isMasked(c);
		bool prev_mask = curr_mask;
		size_t i = beg;
		while (curr_mask == prev_mask && i < seq.size()) {
			prev_mask = curr_mask;
			c = seq[i];
			curr_mask = isMasked(c);
			if (prev_mask != curr_mask){
				return std::pair<size_t,size_t>(beg, i);
			}
			i++;
		}
		return std::pair<size_t,size_t>(beg, i);
	}

	std::vector<std::string> splitOnMask(const std::string& seq) {
		std::vector<std::string> result;
		std::pair<size_t, size_t> reg = Dna::nextToggleMaskedRegion(seq, 0);

		while (reg.first != seq.size()) {
			std::string curr_seq = seq.substr(reg.first, reg.second - reg.first);
			if (curr_seq.size() > 0 && (curr_seq[0] < 97) && curr_seq[0] != 'N') {
				result.push_back(curr_seq);
			}
			reg = Dna::nextToggleMaskedRegion(seq, reg.second);
		}
		
		return result;
	}

	std::vector<SeqInterval> nonmaskedRegions(const std::string& seq) {
		std::vector<SeqInterval> result;
		std::pair<size_t, size_t> reg = Dna::nextToggleMaskedRegion(seq, 0);

		while (reg.first != seq.size()) {
			char reg_char = seq[reg.first];
			if (!isMasked(reg_char)) {
				result.push_back(reg);
			}
			reg = Dna::nextToggleMaskedRegion(seq, reg.second);
		}
		return result;
	}
	std::vector<std::pair<SeqInterval, std::string>> splitOnMaskWithInterval(const std::string& seq) {
		std::vector<std::pair<SeqInterval, std::string>> result;

		std::pair<size_t, size_t> reg = Dna::nextToggleMaskedRegion(seq, 0);

		while (reg.first != seq.size()) {
			std::string curr_seq = seq.substr(reg.first, reg.second - reg.first);
			if (curr_seq.size() > 0 && (curr_seq[0] < 97) && curr_seq[0] != 'N') {
				result.push_back(std::pair<SeqInterval, std::string>(reg, curr_seq));
			}
			reg = Dna::nextToggleMaskedRegion(seq, reg.second);
		}

		return result;
	}

	bool hasMasked(const std::string& seq) {
		for (const char c: seq) {
			if (isMasked(c)) {
				return true;
			}
		}
		return false;
	}

	std::string canonicalKmer(const std::string &kmer) {
		std::string rc = revcom(kmer);
		bool origCanonical = std::lexicographical_compare(
					kmer.cbegin()
					, kmer.cend()
					, rc.cbegin()
					, rc.cend());
		return origCanonical ? kmer : rc;
	}

	std::vector<uint64_t> getHashes(const std::string& seq, size_t kmer_size, size_t hash_n) {
		std::vector<uint64_t> results;
		ntHashIterator itr(seq, hash_n, kmer_size);
		while (itr != itr.end()) {
			for (size_t i = 0; i < hash_n; i++){
				size_t hash_value = (*itr)[i];
				results.push_back(hash_value);
			}
			++itr;
		}
		return results;
	}

	std::vector<uint64_t> getMinimizerHashes(
			const std::string &seq,
			size_t kmer_size,
			size_t hash_n,
			size_t window_size) {

		std::vector<uint64_t> min_hash_per_window{};
		if (window_size > seq.size() || kmer_size > seq.size()) {
			return min_hash_per_window;
		}
		std::vector<uint64_t> hashes = getHashes(seq, kmer_size, hash_n);

		size_t kmers_in_window = window_size - kmer_size + 1;
		size_t windows_in_seq = seq.size() - window_size + 1;

		for (size_t window_idx = 0; window_idx < windows_in_seq; window_idx++) {
			for (size_t offset = 0; offset < hash_n; offset++) {
				uint64_t min_value = hashes[window_idx*hash_n + offset];
				for (size_t i = window_idx; i < window_idx+kmers_in_window; i++){
					min_value = std::min(min_value, hashes[i*hash_n + offset]);
				}
				min_hash_per_window.push_back(min_value);
			}
		}

		return min_hash_per_window;
	}

	void softmaskNotInKmerHashes(std::string& seq, const robin_hood::unordered_set<uint64_t>& kmer_hashes, size_t kmer_size) {
		size_t hash_n = 1;
		for (auto split_rec: Dna::splitOnMaskWithInterval(seq)) {
			if (split_rec.second.size() >= kmer_size) {
				ntHashIterator itr(split_rec.second, hash_n, kmer_size);
				size_t global_beg = split_rec.first.first;
				size_t local_beg = 0;
				while (itr != itr.end()) {
					uint64_t hash_value = (*itr)[0];
					if (!kmer_hashes.count(hash_value)) {
							//std::cerr << "Masking " << rec->seq_id
									//<< '\t' << local_beg  << " .. " << local_beg+kmer_size
									//<< '\t' << rec->seq.substr(global_beg+local_beg, kmer_size)
									//<< '\t' << split_rec.second.substr(local_beg, kmer_size)
									//<< '\t' << hash_value << '\n';
						Dna::softmask(seq, global_beg+local_beg, global_beg+local_beg+kmer_size);
					}
					++itr;
					++local_beg;
				}
			 }
		}
	}
	double shannon(const std::string& seq) {
		const size_t kmer_size = 3;
		if (seq.size() < kmer_size) {
			return 0.0;
		}
		std::map<std::string, int> counts;
		for (size_t i = 0; i < seq.size() - kmer_size; i++) {
			std::string kmer = seq.substr(i, kmer_size);
			counts[kmer]++;
		}
		size_t total_kmer_n = seq.size() - kmer_size + 1;
		double shannon = 0.0;
		for (const auto& pair : counts) {
			double p = static_cast<double>(pair.second) / total_kmer_n;
			shannon += p*std::log(p);
		}
		return -shannon;
	}
	MaskingStats maskingStats(const std::string& seq) {
		MaskingStats result = MaskingStats();
		result.size = seq.size();
		result.softmasked = 0;
		result.hardmasked = 0;
		if (seq.size() == 0) {
			return  result;
		}
		for (const char& c: seq) {
			if (c == 'N' || c == 'n') {
				result.hardmasked++;
			} 
			if (c >= 'a' && c <= 'z' && c != 'n') {
				result.softmasked++;
			}
		}
		return result;
	}
}
