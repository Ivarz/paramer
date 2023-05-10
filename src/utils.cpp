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
		if (buffer == NULL) {
			std::cerr << "Failed to read\n";
			state = ReaderState::ERROR;
		}
		if (retval == NULL) {
			state = ReaderState::DONE;
			std::string line = "";
			return line;
		} else {
			while (!trimNewlineInplace(line)) {
				retval = gzgets(file_handler, buffer, static_cast<int>(buf_size));
				line += std::string(buffer);
			}
			trimNewlineInplace(line);
			last_line = line;
			return line;
		}
	}
	std::optional<std::vector<uint8_t>> Reader::bufferedLoad(uint64_t bytes) {
		int max_bytes = std::numeric_limits<int>::max();
		size_t offset = 0;
		uint64_t loaded_bytes = 0;
		std::vector<uint8_t> bytevec(bytes, 0);

		while (loaded_bytes < bytes) {
			int buffer_size = std::min(static_cast<uint64_t>(max_bytes), bytes - loaded_bytes);
			int gzread_output = gzread(file_handler, reinterpret_cast<char*>(bytevec.data()+offset), buffer_size);
			if (gzread_output < 0) {
				std::cerr << __FUNCTION__ << " Error " << gzerror(file_handler, &gzread_output)  << "\n";
				return  {};
			} else {
				loaded_bytes += buffer_size;
				offset = loaded_bytes;
			}
		}
		std::cout << __FUNCTION__ << " loaded_bytes " << loaded_bytes << '\n';
		return bytevec;

	}
	int Reader::read(void* buff, size_t bytes) {
		return gzread(file_handler, reinterpret_cast<char*>(buff), bytes);
	}
	Writer::Writer(const std::string& fn) : file_name(fn) {
		if (!std::filesystem::exists(file_name)) {
			file_handler = gzopen(file_name.c_str(), "wb");
		} else {
			std::cerr << "Warning overwriting " << file_name << '\n';
			file_handler = gzopen(file_name.c_str(), "wb");
			//file_handler = nullptr;
		}
	}
	Writer::~Writer() {
		if (file_handler) {
			gzclose(file_handler);
		}
	}

	int Writer::write(void* buff, size_t bytes) {
		return gzwrite(file_handler, (char*) buff, bytes);
	}
	int Writer::bufferedWrite(const std::vector<uint8_t>& data) {
		uint64_t written_bytes = 0;
		size_t bytes = data.size();
		int max_bytes = std::numeric_limits<int>::max();
		size_t offset = 0;

		while (written_bytes < bytes) {
			int buffer_size = std::min(static_cast<uint64_t>(max_bytes), bytes - written_bytes);
			int gzwrite_output = gzwrite(file_handler, (char*) &data[offset], buffer_size*sizeof(data.at(offset)));
			if (gzwrite_output < 0) {
				return gzwrite_output;
			} else {
				written_bytes += buffer_size;
				offset = written_bytes;
			}
		}
		return 0;
	}
	int Writer::writeLine(const std::string& str) {
		std::vector<uint8_t> data(str.begin(), str.end());
		data.push_back('\n');
		return bufferedWrite(data);
	}
}
