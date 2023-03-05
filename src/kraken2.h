#ifndef KRAKEN2_H
#define KRAKEN2_H

#include <optional>
#include <vector>
#include <sstream>
#include "utils.h"

namespace Kraken2 {
	using TaxaKmerPair = std::pair<std::string, size_t>;
	class Rec {
		public:
			Rec(bool c
					, const std::string& sid
					, const std::string& tid
					, bool pe
					, size_t r1
					, size_t r2
					, std::vector<TaxaKmerPair> r1_km
					, std::vector<TaxaKmerPair> r2_km
			   );
			void print() const;
			const std::vector<TaxaKmerPair>& getR1Kmers() const { return r1_kmers; }
			std::string seq_id;
		private:
			bool classified;
			std::string taxid;
			bool paired_end;
			size_t r1_size;
			size_t r2_size;
			std::vector<TaxaKmerPair> r1_kmers;
			std::vector<TaxaKmerPair> r2_kmers;
	};

	std::optional<Rec> nextRecord(Gz::Reader& gzr);
}

#endif
