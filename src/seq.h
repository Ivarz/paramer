#ifndef SEQ_H
#define SEQ_H

#include <string>

namespace Dna {
	std::string revcom(const std::string& seq);
	void softmask(std::string& seq, size_t beg, size_t end);
}

#endif
