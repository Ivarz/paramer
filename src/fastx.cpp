#include "fastx.h"
#include <string_view>

namespace Fastx {
	std::optional<FileFormat> inferFileFormat(const std::string& fname) {
		Gz::Reader reader = Gz::Reader(fname);
		std::string curr_line = reader.nextLine();
		if (curr_line[0] == '>') {
			return FileFormat::Fasta;
		} else if (curr_line[0] == '@') {
			return FileFormat::Fastq;
		} else {
			return {};
		}
	}
}
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

	std::optional<Pair> nextRecordPair(Gz::Reader& gzr1, Gz::Reader& gzr2) {
		std::optional<Rec> rec1 = nextRecord(gzr1);
		std::optional<Rec> rec2 = nextRecord(gzr2);
		if (rec1 && rec2) {
			return std::pair<Rec, Rec>(*rec1, *rec2);
		} else {
			return {};
		}
	}

	std::vector<Rec> Rec::splitOnMask() const {
		std::vector<Rec> result;
		std::vector<std::string> split_seqs = Dna::splitOnMask(seq);
		for (size_t i = 0; i < split_seqs.size(); i++) {
			std::string curr_seq_id = seq_id + "_" + std::to_string(i);
			result.emplace_back(curr_seq_id, split_seqs[i], "");
		}
		return result;
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
				Dna::softmask(seq, beg, end);
			}
			kmer_pos += p.second;
		}
		return;
	}

	robin_hood::unordered_set<std::string> loadUnmaskedKmers(const std::string& fname, size_t kmer_size) {
		robin_hood::unordered_set<std::string> result;
		Gz::Reader fh = Gz::Reader(fname);
		std::optional<Fasta::Rec> rec = Fasta::nextRecord(fh);
		while (rec) {
			for (Fasta::Rec seq: rec->splitOnMask()) {
				Dna::addKmers(seq.seq, kmer_size, result);
			}
			std::cerr << "Loaded kmers " << result.size() << '\n';
			rec = Fasta::nextRecord(fh);
		}
		return result;
	}

	robin_hood::unordered_set<uint64_t> loadUnmaskedKmerHashes(const std::string& fname, size_t kmer_size) {
		robin_hood::unordered_set<uint64_t> result;
		size_t hash_n = 1;
		Gz::Reader fh = Gz::Reader(fname);

		std::optional<Fasta::Rec> rec = Fasta::nextRecord(fh);
		while (rec) {
			for (Dna::SeqInterval interval: Dna::nonmaskedRegions(rec->seq)) {
				size_t interval_size = interval.second - interval.first;
				if (interval_size >= kmer_size) {
					ntHashIterator itr(rec->seq.c_str()+interval.first, hash_n, kmer_size);
					size_t hash_i = 0;
					size_t kmer_count = interval_size - kmer_size + 1;
					size_t hashes_in_window = hash_n * kmer_count;
					while (hash_i++ < hashes_in_window) {
						uint64_t hash_value = (*itr)[hash_i % hash_n];
						result.insert(hash_value);
						++itr;
					}
				}
			}
			//for (Fasta::Rec seq: rec->splitOnMask()) {
				//if (seq.size() >= kmer_size) {
					////std::cerr << seq.seq <<'\n';
					//ntHashIterator itr(seq.seq, hash_n, kmer_size);
					//size_t beg = 0;
					//while (itr != itr.end()) {
						//uint64_t hash_value = (*itr)[0];
						//std::string curr_seq = seq.seq.substr(beg, kmer_size);
						//result.insert(hash_value);
						//++itr;
						//++beg;
					//}
				//}
			//}
			rec = Fasta::nextRecord(fh);
		}
		return result;
	}

	void dropKmerHashesFound(const std::string& fname, size_t kmer_size, robin_hood::unordered_set<uint64_t>& kmers) {
		Gz::Reader fh = Gz::Reader(fname);
		std::optional<Fasta::Rec> rec = Fasta::nextRecord(fh);
		size_t hash_n = 1;
		while (rec) {
			for (Fasta::Rec seq: rec->splitOnMask()) {
				std::cerr << "Dropping from " << seq.seq_id << '\n';
				if (seq.size() >= kmer_size) {
					ntHashIterator itr(seq.seq, hash_n, kmer_size);
					while (itr != itr.end()) {
						uint64_t hash_value = (*itr)[0];
						kmers.erase(hash_value);
						++itr;
					}
				}
			}
			rec = Fasta::nextRecord(fh);
		}
	}

	bool idsMatch(const std::string& id1,  const std::string& id2) {
		size_t min_size = std::min(id1.size(), id2.size());
		return id1.substr(0, min_size) == id2.substr(0, min_size);
	}

	void loadSoftmaskAndPrint(const std::string& fasta_fname
			, const std::vector<std::string>& kraken2_fnames
			, const std::vector<std::string>& reference_fnames
			, size_t kmer_size
			) {

		Gz::Reader gzrfa = Gz::Reader(fasta_fname);
		std::vector<Gz::Reader> k2_readers;

		for (const auto &fn : kraken2_fnames) {
			k2_readers.emplace_back(fn);
		}

		std::optional<Fasta::Rec> fa_rec = Fasta::nextRecord(gzrfa);
		std::vector<std::optional<Kraken2::Rec>> k2_recs;

		for (auto &gzkr2 : k2_readers) {
			k2_recs.push_back(Kraken2::nextRecord(gzkr2));
		}

		robin_hood::unordered_set<uint64_t> fa_kmers = {};
		if (reference_fnames.size() > 0) {
			fa_kmers = loadUnmaskedKmerHashes(fasta_fname, kmer_size);
			std::cerr << "kmer hashes loaded: " << fa_kmers.size() << '\n';
			std::cerr << "Dropping hashes found in refs\n";
			for (std::string ref_name : reference_fnames) {
				dropKmerHashesFound(ref_name, kmer_size, fa_kmers);
			}
			std::cerr << "kmer hashes kept: " << fa_kmers.size() << '\n';
			std::cerr << "clean target fasta file\n";
		}

		while (fa_rec) {
			std::cerr << "Processing " << fa_rec->seq_id << '\n';
			for (const auto &k2_rec : k2_recs) {
				if (k2_rec) {
					if (idsMatch(fa_rec->seq_id, k2_rec->seq_id)) {
						std::cerr << "Processing " << fa_rec->seq_id << '\n';
						fa_rec->softmaskWithKraken2(*k2_rec);
						// results.push_back(*fa_rec);
				  } else {
					std::cerr << "Ids do not match\n";
					std::cerr << fa_rec->seq_id << "\t"
							  << k2_rec->seq_id << '\n';
				  }
				}
			}
			if (reference_fnames.size() > 0) {
				Dna::softmaskNotInKmerHashes(fa_rec->seq, fa_kmers, kmer_size);
			}

			fa_rec->print();
			fa_rec = Fasta::nextRecord(gzrfa);
			k2_recs.clear();

			for (auto &gzkr2 : k2_readers) {
			  k2_recs.push_back(Kraken2::nextRecord(gzkr2));
			}
		}
	}
}
