#include "fastx.h"

namespace Fastq {

	Rec::Rec(const std::string& sid, const std::string& sq, const std::string& q) :
		seq_id(sid),
		seq(sq),
		qual(q)
	{}

	void Rec::print() const {
		std::cout << '@' << seq_id << '\n';
		std::cout << seq << '\n';
		std::cout << '+' << '\n';
		std::cout << qual << '\n';
	}

	std::optional<Rec> nextRecord(Gz::Reader& gzr) {
		std::string seq_id = gzr.nextLine();
		std::string seq = gzr.nextLine();
		gzr.nextLine();
		std::string qual = gzr.nextLine();

		if (seq_id.size() > 0) {
			return Rec(seq_id.substr(1), seq, qual);
		} else {
			return {};
		}
	}
}

namespace Fasta {

	Rec::Rec(const std::string& sid, const std::string& sq) :
		seq_id(sid),
		seq(sq)
		{}

		void Rec::print() const {
		std::cout << '>' << seq_id << '\n';
		std::cout << seq << '\n';
	}

	void Rec::softmask(size_t beg, size_t end) {
		Dna::softmask(seq, beg, end);
	}

	std::vector<Rec> Rec::splitOnMask() const {
		std::vector<Rec> result;
		std::vector<std::string> split_seqs = Dna::splitOnMask(seq);
		for (size_t i = 0; i < split_seqs.size(); i++) {
			std::string curr_seq_id = seq_id + "_" + std::to_string(i);
			result.emplace_back(curr_seq_id, split_seqs[i]);
		}
		return result;
	}

	std::optional<Rec> nextRecord(Gz::Reader& gzr) {
		std::string tmp_buffer = gzr.last_line;
		std::string seq_id = tmp_buffer.size() > 0 && tmp_buffer[0] == '>'
			? tmp_buffer : gzr.nextLine();
		trimNewlineInplace(seq_id);
		std::string curr_line = gzr.nextLine();
		std::string seq = "";

		while (curr_line.size() > 0 && curr_line[0] != '>') {
			seq += curr_line;
			curr_line = gzr.nextLine();
		}
		if (seq_id.size() > 0) {
			//seq_id.remove(0);
			return Rec(seq_id.substr(1), seq);
		} else {
			return {};
		}

	}
	void Rec::softmaskWithKraken2(const Kraken2::Rec& k2_rec, size_t kmer_size) {
		auto r1_pairs = k2_rec.getR1Kmers();
		size_t kmer_pos = 0;
		for (const auto& p: r1_pairs) {
			if (p.first != "0"){
				size_t beg = kmer_pos;
				size_t end = kmer_pos + (p.second ? p.second - 1 : p.second) + kmer_size;
				//std::cout << p.first << '\t' << p.second << '\t' << beg << '\t' << end << '\n';
				Dna::softmask(seq, beg, end);
				kmer_pos += p.second;
			}
		}
		return;
	}
}
