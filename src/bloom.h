#ifndef BLOOM_H
#define BLOOM_H
#include <vector>
#include <zlib.h>
#include <fstream>
#include <ntHashIterator.hpp>
#include "seq.h"
#include "fastx.h"

namespace Bloom {
	class Filter {
		public:
			Filter(uint64_t s, uint64_t k , uint64_t h) : filter_size(s), kmer_size(k) , hash_n(h) {
				bytevec = std::vector<uint8_t>(s, 0);
			}
			Filter(const std::string& input_fname);
			void addSeq(const std::string& seq);
			size_t searchSeq(const std::string& seq) const;
			void addFasta(const std::string& fasta_fname, size_t minsize);
			void addFastq(const std::string& fastq_fname, size_t minsize);
			size_t searchFastqPair(const Fastq::Pair& fq_pair) const;
			void write(const std::string& out_fname) const;
			void writeGz(const std::string& out_fname) const;
			size_t size() const { return bytevec.size(); }
		private:
			uint64_t filter_size; // filter size in bytes
			uint64_t kmer_size;
			uint64_t hash_n;
			std::vector<uint8_t> bytevec;
			const size_t ELEMENTS_IN_BUCKET = 8;
			const size_t MAX_BUCKET_VALUE = 255;

	};
}
#endif
