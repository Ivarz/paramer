#ifndef FASTX_H
#define FASTX_H
#include <optional>
#include "utils.h"
#include "seq.h"
#include "kraken2.h"


namespace Fastq {
	struct Rec {
		Rec(const std::string& sid, const std::string& sq, const std::string& q);
		void print() const;
		std::string seq_id;
		std::string seq;
		std::string qual;
	};


	std::optional<Rec> nextRecord(Gz::Reader& gzr);
}

namespace Fasta {
	class Rec {
		public:
			Rec(const std::string& sid, const std::string& sq);

			void print() const;
			void softmask(size_t beg, size_t end);
			void softmaskWithKraken2(const Kraken2::Rec& k2_rec, size_t kmer_size=35);
			size_t size() const { return seq.size(); }
			//const size_t& seq() const { return seq; }
			std::vector<Rec> splitOnMask() const;
			std::string seq_id;
			std::string seq;
	};
	std::optional<Rec> nextRecord(Gz::Reader& gzr);

}

#endif
