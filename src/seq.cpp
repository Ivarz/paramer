#include "seq.h"

namespace Dna {
	std::string revcom(const std::string& seq) {
		std::string result(seq.size(), 'N');
		if (seq.size() == 0) {
			return result;
		}
		for (int i = 0; i < seq.size(); i++) {
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

	std::pair<size_t,size_t> nextToggleMaskedRegion(const std::string& seq, size_t beg) {
		if (beg == seq.size()) {
			return std::pair<size_t,size_t>(beg, beg);
		}
		char c = seq[beg];
		bool curr_mask = c == 'a' || c == 't' || c == 'g' || c == 'c' || c == 'N' || c == 'n';
		bool prev_mask = curr_mask;
		size_t i = beg;
		while (curr_mask == prev_mask && i < seq.size()) {
			prev_mask = curr_mask;
			c = seq[i];
			curr_mask = c == 'a' || c == 't' || c == 'g' || c == 'c' || c == 'N' || c == 'n';
			if (prev_mask != curr_mask){
				return std::pair<size_t,size_t>(beg, i);
			}
			i++;
		}
		return std::pair<size_t,size_t>(beg, i);
	}

	std::vector<std::string> splitOnMask(const std::string& seq) {
		size_t beg = 0;
		size_t end = 0;
		bool masked = false;
		bool prev_masked = masked;
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
}
