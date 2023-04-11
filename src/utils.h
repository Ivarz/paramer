#ifndef UTILS_H
#define UTILS_H
#include <zlib.h>
#include <string>
#include <iostream>
#include <algorithm>
#include <filesystem>

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
		private:
			gzFile file_handler;
	};
}
#endif
