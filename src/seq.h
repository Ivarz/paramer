#ifndef SEQ_H
#define SEQ_H

#include <string>
#include <vector>
#include <set>

namespace Dna {
	std::string revcom(const std::string& seq);
	void softmask(std::string& seq, size_t beg, size_t end);
	void unmask(std::string& seq, size_t beg, size_t end);
	//std::set<std::string> kmers(const std::string& seq, size_t size);
	void addKmers(const std::string& seq, size_t size, std::set<std::string>& kmers);

	std::pair<size_t,size_t> nextToggleMaskedRegion(const std::string& seq, size_t beg);
	std::vector<std::string> splitOnMask(const std::string& seq);
}

#endif
