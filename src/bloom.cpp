#include "bloom.h"
#include "bit_lookup.h"
#include "seq.h"
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <functional>
namespace Bloom {
	std::pair<size_t, uint8_t> index_value(uint64_t hash_value, size_t filter_size) {
		uint64_t bit_idx = hash_value % (filter_size*BITS_IN_BYTE);
		uint64_t byte_idx = bit_idx / BITS_IN_BYTE;
		uint8_t byte_value = 1 << (bit_idx % BITS_IN_BYTE);
		return std::pair<size_t, uint8_t>(byte_idx, byte_value);
	}
	void Filter::addMinimizers(const std::string& seq) {
		for (uint64_t minimizer: Dna::getMinimizerHashes(seq, kmer_size, hash_n, window_size)) {
			std::pair<size_t, uint8_t> idx_value = index_value(minimizer, filter_size);
			bytevec[idx_value.first] |= idx_value.second;
		}
	}
	void Filter::addSeq(const std::string& seq) {
		for (uint64_t minimizer: Dna::getHashes(seq, kmer_size, hash_n)) {
			std::pair<size_t, uint8_t> idx_value = index_value(minimizer, filter_size);
			bytevec[idx_value.first] |= idx_value.second;
		}
		//ntHashIterator itr(seq, hash_n, kmer_size);
		//while (itr != itr.end()) {
			//for (size_t i = 0; i < hash_n; i++){
				//uint64_t hash_value = (*itr)[i];
				//std::pair<size_t, uint8_t> idx_value = index_value(hash_value, filter_size);
				//bytevec[idx_value.first] |= idx_value.second;
			//}
			//++itr;
		//}
	}

	size_t Filter::searchMinimizers(const std::string& seq) const {
		std::vector<uint64_t> minimizers = Dna::getMinimizerHashes(seq, kmer_size, hash_n, window_size);
		size_t kmer_hits = 0;
		for (size_t i = 0; i < minimizers.size(); i += hash_n) {
			size_t hash_hits = 0;
			for (size_t j = i; j < i+hash_n; j++) {
				std::pair<size_t, uint8_t> idx_value = index_value(minimizers[j], filter_size);
				if (bytevec[idx_value.first] & idx_value.second) {
					hash_hits++;
				}
				kmer_hits += (hash_hits/hash_n);
			}
		}
		return kmer_hits;
	}

	size_t Filter::searchSeq(const std::string& seq) const {
		std::vector<uint64_t> hashes = Dna::getHashes(seq, kmer_size, hash_n);
		size_t kmer_hits = 0;
		for (size_t i = 0; i < hashes.size(); i += hash_n) {
			size_t hash_hits = 0;
			for (size_t j = i; j < i+hash_n; j++) {
				std::pair<size_t, uint8_t> idx_value = index_value(hashes[j], filter_size);
				if (bytevec[idx_value.first] & idx_value.second) {
					hash_hits++;
				}
			kmer_hits += (hash_hits/hash_n);
		}
		
		//ntHashIterator itr(seq, hash_n, kmer_size);
		//while (itr != itr.end()) {
			//size_t hash_hits = 0;
			//for (size_t i = 0; i < hash_n; i++){
				//uint64_t hash_value = (*itr)[i];
				//std::pair<size_t, uint8_t> idx_value = index_value(hash_value, filter_size);
				//if (bytevec[idx_value.first] & idx_value.second) {
					//hash_hits++;
				//}
			//}
			//if (hash_hits == hash_n) {
				//kmer_hits++;
			//}
			//++itr;
		}
		return kmer_hits;
	}

	size_t Filter::searchFastqPair(const Fastq::Pair& fq_pair) const {
		size_t kmer_hits = searchSeq(fq_pair.first.seq);
		kmer_hits += searchSeq(fq_pair.second.seq);
		return kmer_hits;
	}

	void Filter::addFasta(const std::string& fasta_fname, size_t minsize) {
		Gz::Reader gzrfa(fasta_fname);
		std::optional<Fasta::Rec> fa_rec = Fasta::nextRecord(gzrfa);
		while (fa_rec) {
			if (fa_rec) {
				for (const auto& rec: fa_rec->splitOnMask()) {
					if (rec.size() >= minsize) {
						addSeq(rec.seq);
					}
				}
			}
			fa_rec = Fasta::nextRecord(gzrfa);
		}
	}

	void Filter::addFastq(const std::string& fastq_fname, size_t minsize) {
		Gz::Reader gzrfq(fastq_fname);
		std::optional<Fastq::Rec> fq_rec = Fastq::nextRecord(gzrfq);
		size_t counter = 0;
		while (fq_rec) {
			if (fq_rec) {
				addSeq(fq_rec->seq);
				//for (const auto& rec: fq_rec->splitOnMask()) {
					//if (rec.size() >= minsize) {
						//addSeq(rec.seq);
					//}
				//}
			}
			fq_rec = Fastq::nextRecord(gzrfq);
			counter++;
			if (!(counter % 1000000)) {
				std::cerr << counter << '\n';
				std::cerr << fq_rec->seq << '\n';
			}
		}
	}


	void Filter::writeRaw(const std::string& out_fname) const {
		std::ofstream outfh(out_fname, std::ios::out | std::ios::binary);
		uint64_t fsize = static_cast<uint64_t>(bytevec.size());
		uint64_t k = static_cast<uint64_t>(kmer_size);
		uint64_t w = static_cast<uint64_t>(window_size);
		uint64_t h = static_cast<uint64_t>(hash_n);
		outfh.write((char*) &fsize, sizeof(fsize));
		outfh.write((char*) &k, sizeof(k));
		outfh.write((char*) &w, sizeof(w));
		outfh.write((char*) &h, sizeof(h));
		//outfh.write((char*) &out_fname[0], out_fname.size()*sizeof(char));
		outfh.write((char*) &bytevec[0], filter_size*sizeof(bytevec.at(0)));
		outfh.close();
	}

	int Filter::write(const std::string& out_fname) const {
		Gz::Writer gzwriter(out_fname);
		//gzFile fp = gzopen(out_fname.c_str(),"wb");

		uint64_t fsize = static_cast<uint64_t>(bytevec.size());
		uint64_t k = static_cast<uint64_t>(kmer_size);
		uint64_t w = static_cast<uint64_t>(window_size);
		uint64_t h = static_cast<uint64_t>(hash_n);
		gzwriter.write(&fsize, sizeof(fsize));
		gzwriter.write(&k, sizeof(k));
		gzwriter.write(&k, sizeof(w));
		gzwriter.write(&h, sizeof(h));
		return gzwriter.bufferedWrite(bytevec);
	}

	std::optional<Filter> Filter::loadRaw(const std::string& in_fname) {
		uint64_t filter_size = 0;
		uint64_t kmer_size = 0;
		uint64_t window_size = 0;
		uint64_t hash_n = 0;

		std::ifstream infh(in_fname, std::ios::out | std::ios::binary);
		infh.read(reinterpret_cast<char*>(&filter_size), sizeof(uint64_t));
		infh.read(reinterpret_cast<char*>(&kmer_size), sizeof(uint64_t));
		infh.read(reinterpret_cast<char*>(&window_size), sizeof(uint64_t));
		infh.read(reinterpret_cast<char*>(&hash_n), sizeof(uint64_t));

		// Read the remaining data into a vector<uint8_t>
		Filter result = Filter(filter_size, kmer_size, window_size, hash_n);
		infh.read(reinterpret_cast<char*>(result.bytevec.data()), filter_size);

		// Close the file
		infh.close();
		return result;
	}

	std::optional<Filter> Filter::load(const std::string& in_fname) {

		uint64_t filter_size = 0;
		uint64_t kmer_size = 0;
		uint64_t window_size = 0;
		uint64_t hash_n = 0;

		Gz::Reader gz_reader(in_fname);

		gz_reader.read(&filter_size, sizeof(uint64_t));
		gz_reader.read(&kmer_size, sizeof(uint64_t));
		gz_reader.read(&window_size, sizeof(uint64_t));
		gz_reader.read(&hash_n, sizeof(uint64_t));

		// Read the remaining data into a vector<uint8_t>
		auto bufload = gz_reader.bufferedLoad(filter_size);
		if (bufload) {
			Filter result = Filter(0, 0, 0, 0);
			result.filter_size = filter_size;
			result.kmer_size = kmer_size;
			result.hash_n = hash_n;
			result.bytevec = std::move(*bufload);
			return result;
			//return std::move(result);
		} else {
			return {};
		}
	}
	size_t Filter::setBitsCount() const {
		size_t count = 0;
		for (size_t i=0; i < size(); i++) {
			count += BitLookup::BIT_COUNT[bytevec[i]];
		}
		return count;
	}
	double Filter::falsePostiveRate() const {
		double k = static_cast<double>(hash_n);
		double m = static_cast<double>(bytevec.size());
		double set_bits = static_cast<double>(this->setBitsCount());
		double n_star = - (m/k)*logf(1.0-(set_bits/m));
		double res = pow(1 - exp((-k*n_star)/m), k);
		return res;
	}
	void Filter::dfs(std::string current_seq, robin_hood::unordered_set<uint64_t> seen_kmer_hashes,
                 std::vector<std::string>& candidate_seqs,
                 const std::function<std::string(std::string)>& extract_kmer,
                 const std::function<std::string(std::string, char)>& next_seq
				 ) const {
		std::string current_kmer = extract_kmer(current_seq);
		Dna::addKmerHashes(current_kmer, kmer_size, 1, seen_kmer_hashes);
		bool finished_path = true;
		for (char c : std::string("ATGC")) {
			std::string potential_neighbour = extract_kmer(next_seq(current_kmer, c));
			uint64_t neighbour_hash = Dna::getHashes(potential_neighbour, kmer_size, 1)[0];
			if (searchSeq(potential_neighbour) && !seen_kmer_hashes.count(neighbour_hash)) {
				finished_path = false;
				dfs(next_seq(current_seq, c), seen_kmer_hashes, candidate_seqs, extract_kmer, next_seq);
			}
		}
		if (finished_path) {
			//std::cout << "Candidate\t" << current_seq << '\n';
			candidate_seqs.push_back(current_seq);
		}
		//std::cout << "Backtrack\t" << current_seq << '\n';
	}

	void Filter::dfs5prime(const std::string& current_seq,
			robin_hood::unordered_set<uint64_t> seen_kmer_hashes,
			std::vector<std::string>& candidate_seqs
			) const {

		auto extract_kmer = [this](std::string s){ return s.substr(s.size()-kmer_size, kmer_size); };
		auto next_seq = [](std::string str, char c) { return str + c; };
		dfs(current_seq,
				seen_kmer_hashes,
				candidate_seqs,
				extract_kmer,
				next_seq
				);
	}
	void Filter::dfs3prime(const std::string& current_seq,
			robin_hood::unordered_set<uint64_t> seen_kmer_hashes,
			std::vector<std::string>& candidate_seqs
			) const {

		auto extract_kmer = [this](std::string s){ return s.substr(0, kmer_size); };
		auto next_seq = [](std::string str, char c) { return c + str; };
		dfs(current_seq,
				seen_kmer_hashes,
				candidate_seqs,
				extract_kmer,
				next_seq
				);
	}
	std::vector<std::string> Filter::extendSeq(const std::string &seq) const {
		robin_hood::unordered_set<uint64_t> seen_kmers_5p;
		size_t hash_n_for_dfs = 1;
		Dna::addKmerHashes(seq, kmer_size, hash_n_for_dfs, seen_kmers_5p);
		robin_hood::unordered_set<uint64_t> seen_kmers_3p = seen_kmers_5p;
		std::vector<std::string> candidate_seqs_5p;
		std::vector<std::string> candidate_seqs_3p;
		std::vector<std::string> extended_seqs;

		dfs5prime(seq, seen_kmers_5p, candidate_seqs_5p);
		dfs3prime(seq, seen_kmers_3p, candidate_seqs_3p);
		for (const auto& seq5: candidate_seqs_5p) {
			for (const auto& seq3: candidate_seqs_3p) {
				std::string seq5_extension = seq5.substr(seq.size(), seq5.size()-seq.size());
				std::string seq3_extension = seq3.substr(0, seq3.size()-seq.size());
				extended_seqs.push_back(seq3_extension + seq + seq5_extension);
			}
		}

		return extended_seqs;
	}

}
