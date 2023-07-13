#include "bloom.h"
#include "bit_lookup.h"
#include "seq.h"
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <functional>
#include <queue>
#include <map>
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
	}

	size_t Filter::searchMinimizers(const std::string& seq) {
		if (seq.size() < window_size || seq.find('N') != std::string::npos || seq.find('n') != std::string::npos) {
			return 0;
		}
		std::vector<uint64_t> minimizers = Dna::getMinimizerHashes(seq, kmer_size, hash_n, window_size);
		size_t kmer_hits = 0;
		for (size_t i = 0; i < minimizers.size(); i += hash_n) {
			size_t hash_hits = 0;
			for (size_t j = i; j < i+hash_n; j++) {
				std::pair<size_t, uint8_t> idx_value = index_value(minimizers[j], filter_size);
				//if (bytevec[idx_value.first] & idx_value.second) {
				if (getByteVecVal(idx_value.first) & idx_value.second) {
					hash_hits++;
				}
				kmer_hits += (hash_hits/hash_n);
			}
		}
		return kmer_hits;
	}

	size_t Filter::searchSeq(const std::string& seq) {
		if (seq.size() < kmer_size || seq.find('N') != std::string::npos || seq.find('n') != std::string::npos) {
			return 0;
		}
		std::vector<uint64_t> hashes = Dna::getHashes(seq, kmer_size, hash_n);
		size_t kmer_hits = 0;
		for (size_t i = 0; i < hashes.size(); i += hash_n) {
			size_t hash_hits = 0;
			for (size_t j = i; j < i+hash_n; j++) {
				std::pair<size_t, uint8_t> idx_value = index_value(hashes[j], filter_size);
				if (getByteVecVal(idx_value.first) & idx_value.second) {
					hash_hits++;
				}
			}
			kmer_hits += (hash_hits/hash_n);
		}
		return kmer_hits;
	}

	uint8_t Filter::seekAt(size_t idx) {
		uint8_t value = 0;
		const size_t OFFSET = 34; //header subject to change
		filter_fptr.seekg(OFFSET+idx, filter_fptr.beg);
		filter_fptr.read(reinterpret_cast<char*>(&value), sizeof(value));

		return value;
	}

	size_t Filter::seekSeq(const std::string& seq) {
		std::vector<uint64_t> hashes = Dna::getHashes(seq, kmer_size, hash_n);
		size_t kmer_hits = 0;
		for (size_t i = 0; i < hashes.size(); i += hash_n) {
			size_t hash_hits = 0;
			for (size_t j = i; j < i+hash_n; j++) {
				std::pair<size_t, uint8_t> idx_value = index_value(hashes[j], filter_size);
				if (seekAt(idx_value.first) & idx_value.second) {
					hash_hits++;
				}
			}
			kmer_hits += (hash_hits/hash_n);
		}
		return kmer_hits;
	}

	size_t Filter::searchFastqPair(const Fastq::Pair& fq_pair) {

		size_t (Filter::*searcher)(const std::string&);
		searcher = 
			window_size > kmer_size
			? &Filter::searchMinimizers
			: &Filter::searchSeq;
		size_t kmer_hits = (this->*searcher)(fq_pair.first.seq);
		kmer_hits += (this->*searcher)(fq_pair.second.seq);
		//size_t kmer_hits = searchSeq(fq_pair.first.seq);
		//kmer_hits += searchSeq(fq_pair.second.seq);
		return kmer_hits;
	}

	void Filter::addFasta(const std::string& fasta_fname, size_t minsize) {
		Gz::Reader gzrfa(fasta_fname);
		std::optional<Fasta::Rec> fa_rec = Fasta::nextRecord(gzrfa);
		void (Filter::*adder)(const std::string&);
		adder = window_size > kmer_size ? &Filter::addMinimizers : &Filter::addSeq;
		while (fa_rec) {
			if (fa_rec) {
				for (const auto& rec: fa_rec->splitOnMask()) {
					if (rec.size() >= minsize) {
						//addSeq(rec.seq);
						(this->*adder)(rec.seq);
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
		void (Filter::*adder)(const std::string&);
		adder = window_size > kmer_size ? &Filter::addMinimizers : &Filter::addSeq;
		while (fq_rec) {
			if (fq_rec) {
				//addSeq(fq_rec->seq);
				(this->*adder)(fq_rec->seq);
			}
			fq_rec = Fastq::nextRecord(gzrfq);
			counter++;
			if (!(counter % 1000000)) {
				std::cerr << counter << '\n';
				std::cerr << fq_rec->seq << '\n';
			}
		}
	}


	int Filter::writeRaw(const std::string& out_fname) const {
		std::ofstream outfh(out_fname, std::ios::out | std::ios::binary);
		uint8_t magic_byte = 0;
		uint64_t fsize = static_cast<uint64_t>(bytevec.size());
		uint64_t k = static_cast<uint64_t>(kmer_size);
		uint64_t w = static_cast<uint64_t>(window_size);
		uint64_t h = static_cast<uint64_t>(hash_n);
		outfh.write((char*) &magic_byte, sizeof(magic_byte));
		outfh.write((char*) &magic_byte, sizeof(magic_byte));
		outfh.write((char*) &fsize, sizeof(fsize));
		outfh.write((char*) &k, sizeof(k));
		outfh.write((char*) &w, sizeof(w));
		outfh.write((char*) &h, sizeof(h));
		//outfh.write((char*) &out_fname[0], out_fname.size()*sizeof(char));
		outfh.write((char*) &bytevec[0], filter_size*sizeof(bytevec.at(0)));
		outfh.close();
		return 0;
	}

	int Filter::writeGz(const std::string& out_fname) const {
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

	int Filter::write(const std::string& out_fname, Compression cmpr) const {
		if (cmpr == Compression::GZ) {
			return writeGz(out_fname);
		} else {
			return writeRaw(out_fname);
		}
	}

	std::optional<Filter> Filter::loadRaw(const std::string& in_fname) {
		uint8_t magic_byte = 0;
		uint64_t filter_size = 0;
		uint64_t kmer_size = 0;
		uint64_t window_size = 0;
		uint64_t hash_n = 0;

		std::ifstream infh(in_fname, std::ios::out | std::ios::binary);
		infh.read(reinterpret_cast<char*>(&magic_byte), sizeof(magic_byte));
		infh.read(reinterpret_cast<char*>(&magic_byte), sizeof(magic_byte));
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

	std::optional<Filter> Filter::loadPointer(const std::string& in_fname) {
		uint8_t magic_byte = 0;
		uint64_t filter_size = 0;
		uint64_t kmer_size = 0;
		uint64_t window_size = 0;
		uint64_t hash_n = 0;

		std::ifstream infh(in_fname, std::ios::out | std::ios::binary);
		infh.read(reinterpret_cast<char*>(&magic_byte), sizeof(magic_byte));
		infh.read(reinterpret_cast<char*>(&magic_byte), sizeof(magic_byte));
		infh.read(reinterpret_cast<char*>(&filter_size), sizeof(uint64_t));
		infh.read(reinterpret_cast<char*>(&kmer_size), sizeof(uint64_t));
		infh.read(reinterpret_cast<char*>(&window_size), sizeof(uint64_t));
		infh.read(reinterpret_cast<char*>(&hash_n), sizeof(uint64_t));

		// Read the remaining data into a vector<uint8_t>
		Filter result = Filter(0, 0, 0, 0);
		result.filter_size = filter_size;
		result.kmer_size = kmer_size;
		result.window_size = window_size;
		result.hash_n = hash_n;
		result.bytevec = {};
		result.filter_fptr.open(in_fname, std::ios::out | std::ios::binary);

		// Close the file
		infh.close();
		return result;
	}

	void Filter::closePointer(){
		if (filter_fptr.is_open()) {
			filter_fptr.close();
		}
	}

	std::optional<Filter> Filter::loadGz(const std::string& in_fname) {

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

  
	std::optional<Filter> Filter::load(const std::string& in_fname) {
		Bloom::Compression cmpr = Filter::inferCompression(in_fname);
		if (cmpr == Bloom::Compression::GZ) {
			return Filter::loadGz(in_fname);
		} else {
			return Filter::loadRaw(in_fname);
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

	void Filter::bfs(std::string src, robin_hood::unordered_set<uint64_t> seen_kmer_hashes,
                 std::vector<std::string>& candidate_seqs,
                 const std::function<std::string(std::string)>& extract_kmer,
                 const std::function<std::string(std::string, char)>& next_seq,
                 const std::function<std::string(std::string, std::string)>& add_parent_to_path,
                 const std::function<std::string(std::string, uint64_t)>& clip
				 ) {
		std::string current_kmer = extract_kmer(src);
		Dna::addKmerHashes(current_kmer, kmer_size, 1, seen_kmer_hashes);
		bool finished_path = true;
		std::queue<std::string> q{};
		std::map<std::string, std::string> parent{};
		std::map<std::string, int> pathlen{};
		q.push(current_kmer);
		parent[current_kmer] = "";
		pathlen[current_kmer] = 0;

		while (q.size() > 0) {
			std::string curr_seq = q.front();
			//std::cerr << curr_seq << '\t' << pathlen[curr_seq] << "\tseen kmers: " << seen_kmer_hashes.size() << '\n';
			q.pop();
			bool is_leaf = true;
			for (char c : std::string("ATGC")) {
				std::string potential_neighbour = extract_kmer(next_seq(curr_seq, c));
				uint64_t neighbour_hash = Dna::getHashes(potential_neighbour, kmer_size, 1)[0];
				if (searchSeq(potential_neighbour) && !seen_kmer_hashes.count(neighbour_hash)) {
					Dna::addKmerHashes(potential_neighbour, kmer_size, 1, seen_kmer_hashes);
					parent[potential_neighbour] = curr_seq;
					pathlen[potential_neighbour] = pathlen[curr_seq]+1;
					q.push(potential_neighbour);
					is_leaf = false;
				}
			}
			if (is_leaf) {
				std::string path = curr_seq;
				std::string curr_parent = parent[path];
				while (curr_parent != "") {
					path = add_parent_to_path(path, curr_parent);
					curr_parent = parent[curr_parent];
				}
				candidate_seqs.push_back(clip(path, kmer_size));
			}
		}
	}

	void Filter::bfs5prime(const std::string& current_seq,
			robin_hood::unordered_set<uint64_t> seen_kmer_hashes,
			std::vector<std::string>& candidate_seqs
			) {

		auto extract_kmer = [this](std::string s){ return s.substr(s.size()-kmer_size, kmer_size); };
		auto next_seq = [](std::string str, char c) { return str + c; };
		auto add_parent_to_path = [](std::string path, std::string parent) { return parent[0] + path; };
		auto clip = [](std::string path, uint64_t ksize) { return path.substr(ksize, path.size()-ksize); };
		bfs(current_seq,
				seen_kmer_hashes,
				candidate_seqs,
				extract_kmer,
				next_seq,
				add_parent_to_path,
				clip
				);
	}

	void Filter::bfs3prime(const std::string& current_seq,
			robin_hood::unordered_set<uint64_t> seen_kmer_hashes,
			std::vector<std::string>& candidate_seqs
			) {

		auto extract_kmer = [this](std::string s){ return s.substr(0, kmer_size); };
		auto next_seq = [](std::string str, char c) { return c + str; };
		auto add_parent_to_path = [](std::string path, std::string parent) { return  path + parent[parent.size()-1]; };
		auto clip = [](std::string path, uint64_t ksize) { return path.substr(0, path.size()-ksize); };
		bfs(current_seq,
				seen_kmer_hashes,
				candidate_seqs,
				extract_kmer,
				next_seq,
				add_parent_to_path,
				clip
				);
	}

	void Filter::dfs(std::string current_seq, robin_hood::unordered_set<uint64_t> seen_kmer_hashes,
                 std::vector<std::string>& candidate_seqs,
                 const std::function<std::string(std::string)>& extract_kmer,
                 const std::function<std::string(std::string, char)>& next_seq
				 ) {
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
			) {

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
			) {

		auto extract_kmer = [this](std::string s){ return s.substr(0, kmer_size); };
		auto next_seq = [](std::string str, char c) { return c + str; };
		dfs(current_seq,
				seen_kmer_hashes,
				candidate_seqs,
				extract_kmer,
				next_seq
				);
	}
	std::vector<std::string> Filter::extendSeq(const std::string &seq) {
		robin_hood::unordered_set<uint64_t> seen_kmers_5p;
		size_t hash_n_for_dfs = 1;
		Dna::addKmerHashes(seq, kmer_size, hash_n_for_dfs, seen_kmers_5p);
		robin_hood::unordered_set<uint64_t> seen_kmers_3p = seen_kmers_5p;
		std::vector<std::string> candidate_seqs_5p;
		std::vector<std::string> candidate_seqs_3p;
		std::vector<std::string> extended_seqs;

		//dfs5prime(seq, seen_kmers_5p, candidate_seqs_5p);
		//dfs3prime(seq, seen_kmers_3p, candidate_seqs_3p);
		//std::cerr << "bfs5prime\n";
		bfs5prime(seq, seen_kmers_5p, candidate_seqs_5p);
		bfs3prime(seq, seen_kmers_3p, candidate_seqs_3p);
		//std::cerr << "candidate sizes: " << candidate_seqs_5p.size() << '\t' << candidate_seqs_3p.size() << '\n';
		for (const auto& seq5: candidate_seqs_5p) {
			for (const auto& seq3: candidate_seqs_3p) {
				//std::cerr << "seq3 " << seq3 << '\n';
				//std::cerr << "seq5 " << seq5 << '\n';
				//std::string seq5_extension = seq5.substr(seq.size(), seq5.size()-seq.size());
				//std::string seq3_extension = seq3.substr(0, seq3.size()-seq.size());
				//extended_seqs.push_back(seq3_extension + seq + seq5_extension);
				extended_seqs.push_back(seq3 + seq + seq5);
			}
		}
		for (const auto& s: extended_seqs) {
			//std::cerr << "ext seq " << s << '\n';
		}

		return extended_seqs;
	}

	std::vector<std::string> Filter::extendSeqPair(const std::string &seq1,const std::string &seq2 ) {
		std::vector<std::string> result{};
		std::vector<std::string> candidates_mate1 = extendSeq(seq1);
		std::vector<std::string> candidates_mate2 = extendSeq(seq2);
		const std::string seq1_rc = Dna::revcom(seq1);
		const std::string seq2_rc = Dna::revcom(seq2);
		for (auto c1: candidates_mate1) {
			if (c1.find(seq2) != std::string::npos ||c1.find(seq2_rc) != std::string::npos) {
				result.push_back(c1);
			}
		}
		for (auto c2: candidates_mate1) {
			if (c2.find(seq1) != std::string::npos ||c2.find(seq1_rc) != std::string::npos) {
				result.push_back(c2);
			}
		}
		return result;
	}

	Bloom::Compression Filter::inferCompression(const std::string &in_fname) {
		uint8_t magic_byte1 = 0;
		uint8_t magic_byte2 = 0;

		std::ifstream infh(in_fname, std::ios::out | std::ios::binary);
		infh.read(reinterpret_cast<char*>(&magic_byte1), sizeof(magic_byte1));
		infh.read(reinterpret_cast<char*>(&magic_byte2), sizeof(magic_byte2));
		if (magic_byte1 == 0x1f && magic_byte2 == 0x8b) {
			return Bloom::Compression::GZ;
		} else {
			return Bloom::Compression::RAW;
		}
	}

}
