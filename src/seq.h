#ifndef SEQ_H
#define SEQ_H

#include <string>
#include <vector>

namespace Dna {
	std::string revcom(const std::string& seq);
	void softmask(std::string& seq, size_t beg, size_t end);
	void unmask(std::string& seq, size_t beg, size_t end);
	std::pair<size_t,size_t> nextToggleMaskedRegion(const std::string& seq, size_t beg);
	std::vector<std::string> splitOnMask(const std::string& seq);
}

#endif
