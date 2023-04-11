#ifndef UTILS_H
#define UTILS_H
#include <zlib.h>
#include <string>
#include <iostream>
#include <algorithm>
#include <filesystem>
#include <vector>
#include <optional>

bool trimNewlineInplace(std::string& str);

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

			const size_t buf_size = 4096;
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
			~Writer();
			void bufferedWrite(std::vector<uint8_t>& data);
		private:
			std::string file_name;
			gzFile file_handler;
	};
}
#endif
