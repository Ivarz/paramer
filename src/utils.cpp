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
		if (std::filesystem::exists(file_name)) {
			file_handler = gzopen(file_name.c_str(), "rb");
			state = ReaderState::OK;
		} else {
			std::cerr << "File does not exist " << file_name << '\n';
			file_handler = nullptr; //gzopen(file_name.c_str(), "rb");
			state = ReaderState::FILENOTFOUND;
		}
	}

	Reader::Reader(const Reader& gzr) noexcept {
		//std::cerr << "gz cptor\n";
		last_line = gzr.last_line;
		file_name = gzr.file_name;
		file_handler = gzopen(file_name.c_str(), "rb");
		buffer = new char[buf_size];
		state = gzr.state;
	}

	Reader::Reader(Reader&& other) noexcept {
		//std::cerr << "gz mvtor\n";
		last_line = std::move(other.last_line);
		file_name = std::move(other.file_name);
		file_handler = gzopen(file_name.c_str(), "rb");
		buffer = new char[buf_size];
		state = std::move(other.state);
	}


	Reader::~Reader() {
		//std::cerr << "gz dtor\n";
		gzclose(file_handler);
		delete[] buffer;
	}
	std::string Reader::nextLine() {
		if (state == ReaderState::FILENOTFOUND) {
			return "";
		}
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
	//Writer::Writer(const std::string& fn) {
		//uint64_t fsize = static_cast<uint64_t>(bytevec.size());
		//uint64_t k = static_cast<uint64_t>(kmer_size);
		//uint64_t h = static_cast<uint64_t>(hash_n);
		//gzwrite(fp, (char*) &fsize, sizeof(fsize));
		//gzwrite(fp, (char*) &k, sizeof(k));
		//gzwrite(fp, (char*) &h, sizeof(h));
	//}
}
