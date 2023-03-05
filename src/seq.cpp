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

}
