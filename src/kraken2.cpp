#include "kraken2.h"

namespace Kraken2 {
	Rec::Rec(bool c
			, const std::string& sid
			, const std::string& tid
			, bool pe
			, size_t r1
			, size_t r2
			, std::vector<TaxaKmerPair> r1_km
			, std::vector<TaxaKmerPair> r2_km
	   ) :
		classified(c)
		, seq_id(sid)
		, taxid(tid)
		, paired_end(pe)
		, r1_size(r1)
		, r2_size(r2)
		, r1_kmers(r1_km)
		, r2_kmers(r2_km)
	{}

	void Rec::print() const {
		std::cout << classified
			<< '\t' << seq_id
			<< '\t' << taxid
			<< '\t' << paired_end
			<< '\t' << r1_size
			<< '\t' << r2_size;
		for (const auto& p: r1_kmers) {
			std::cout << '\t' << p.first << ':' << p.second;
		}
		std::cout << '\n';
	}


	std::optional<Rec> nextRecord(Gz::Reader& gzr) {
		//Assume single end for current purposes
		const std::string& line = gzr.nextLine();
		if (line.size() > 0) {
			std::istringstream iss(line);
			std::string str, clsf_str, seq_id, taxid, seq_len;

			std::vector<TaxaKmerPair> r1_kmers;

			iss >> clsf_str >> seq_id >> taxid >> seq_len;
			while (iss >> str) {
				size_t pos = str.find(':');
				std::string taxid = str.substr(0, pos);
				size_t kmers = static_cast<size_t>(std::stoi(str.substr(pos+1)));
				TaxaKmerPair kmer_pair(taxid, kmers);
				r1_kmers.push_back(kmer_pair);
			}

			bool clsf = clsf_str == "C";
			return Rec(clsf
					, seq_id
					, taxid
					, false
					, static_cast<size_t>(std::stoi(seq_len))
					, 0
					, r1_kmers
					, {});
		} else {
			return {};
		}

	}
	//Rec nextRecord(Gz::Reader& gzr) {
}

