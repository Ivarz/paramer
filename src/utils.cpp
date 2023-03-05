#include "utils.h"

bool trimNewlineInplace(std::string& str) {
	size_t input_size = str.size();
	str.erase(
			std::remove(str.begin(), str.end(), '\n')
			, str.end());
	str.erase(
			std::remove(str.begin(), str.end(), '\r')
			, str.end());
	size_t output_size = str.size();
	return output_size != input_size;
}

namespace Gz {
	Reader::Reader(const std::string& fn) : file_name(fn) {
		buffer = new char[buf_size];
		file_handler = gzopen(file_name.c_str(), "rb");
	}

	Reader::~Reader() {
		gzclose(file_handler);
		delete[] buffer;
	}
	std::string Reader::nextLine() {
		auto retval = gzgets(file_handler, buffer, static_cast<int>(buf_size));
		std::string line(buffer);
		while (!trimNewlineInplace(line)) {
			retval = gzgets(file_handler, buffer, static_cast<int>(buf_size));
			line += std::string(buffer);
		}
		if (buffer == NULL) {
			std::cerr << "Failed to read\n";
			state = ReaderState::ERROR;
		}
		if (retval == NULL) {
			state = ReaderState::DONE;
			std::string line = "";
			return line;
		} else {
			trimNewlineInplace(line);
			last_line = line;
			return line;
		}
	}


}
