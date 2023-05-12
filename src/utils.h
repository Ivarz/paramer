#ifndef UTILS_H
#define UTILS_H
#include <zlib.h>
#include <string>
#include <iostream>
#include <algorithm>
#include <filesystem>
#include <vector>
#include <optional>

namespace Utils {
	bool trimNewlineInplace(std::string& str);
	size_t dataSizeToBytes(std::string str);
}

namespace Gz {
	enum class ReaderState {
		OK,
		DONE,
		ERROR,
		FILENOTFOUND,
	};

	class Reader {
		public:

			explicit Reader(const std::string& fn);
			Reader(Reader&& other) noexcept;
			Reader(const Reader&) noexcept;
			~Reader();

			std::string nextLine();
			std::optional<std::vector<uint8_t>> bufferedLoad(uint64_t bytes);
			int read(void* buff, size_t bytes);
			std::string last_line = "";

		private:
			Reader& operator=(const Reader&);

			const size_t buf_size = 8192;
			std::string file_name;
			gzFile file_handler;
			char* buffer;
			ReaderState state = ReaderState::OK;

	};
	class Writer {
		public:
			explicit Writer(const std::string& fn);
			Writer(Reader&& other) noexcept;
			Writer(const Reader&) noexcept;
			int write(void* buff, size_t bytes);
			~Writer();
			int bufferedWrite(const std::vector<uint8_t>& data);
			int writeLine(const std::string& str);
		private:
			std::string file_name;
			gzFile file_handler;
	};
}
#endif
