#include "bloom.h"
#include "bit_lookup.h"

namespace Bloom {
	std::pair<size_t, uint8_t> index_value(uint64_t hash_value, size_t filter_size) {
		uint64_t bit_idx = hash_value % (filter_size*BITS_IN_BYTE);
		uint64_t byte_idx = bit_idx / BITS_IN_BYTE;
		uint8_t byte_value = 1 << (bit_idx % BITS_IN_BYTE);
		return std::pair<size_t, uint8_t>(byte_idx, byte_value);
	}
	void Filter::addSeq(const std::string& seq) {
		ntHashIterator itr(seq, hash_n, kmer_size);
		while (itr != itr.end()) {
			for (size_t i = 0; i < hash_n; i++){
				uint64_t hash_value = (*itr)[i];
				std::pair<size_t, uint8_t> idx_value = index_value(hash_value, filter_size);
				bytevec[idx_value.first] |= idx_value.second;
			}
			++itr;
		}
	}

	size_t Filter::searchSeq(const std::string& seq) const {
		ntHashIterator itr(seq, hash_n, kmer_size);
		size_t kmer_hits = 0;
		
		while (itr != itr.end()) {
			size_t hash_hits = 0;
			for (size_t i = 0; i < hash_n; i++){
				uint64_t hash_value = (*itr)[i];
				std::pair<size_t, uint8_t> idx_value = index_value(hash_value, filter_size);
				if (bytevec[idx_value.first] & idx_value.second) {
					hash_hits++;
				}
			}
			if (hash_hits == hash_n) {
				kmer_hits++;
			}
			++itr;
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


	void Filter::write(const std::string& out_fname) const {
		std::ofstream outfh(out_fname, std::ios::out | std::ios::binary);
		uint64_t fsize = static_cast<uint64_t>(bytevec.size());
		uint64_t k = static_cast<uint64_t>(kmer_size);
		uint64_t h = static_cast<uint64_t>(hash_n);
		outfh.write((char*) &fsize, sizeof(fsize));
		outfh.write((char*) &k, sizeof(k));
		outfh.write((char*) &h, sizeof(h));
		//outfh.write((char*) &out_fname[0], out_fname.size()*sizeof(char));
		outfh.write((char*) &bytevec[0], filter_size*sizeof(bytevec.at(0)));
		outfh.close();
	}

	void Filter::writeGz(const std::string& out_fname) const {
		gzFile fp = gzopen(out_fname.c_str(),"wb");

		uint64_t fsize = static_cast<uint64_t>(bytevec.size());
		uint64_t k = static_cast<uint64_t>(kmer_size);
		uint64_t h = static_cast<uint64_t>(hash_n);
		gzwrite(fp, (char*) &fsize, sizeof(fsize));
		gzwrite(fp, (char*) &k, sizeof(k));
		gzwrite(fp, (char*) &h, sizeof(h));

		//outfh.write((char*) &out_fname[0], out_fname.size()*sizeof(char));
		gzwrite(fp, (char*) &bytevec[0], filter_size*sizeof(bytevec.at(0)));

		gzclose(fp);
	}

	//Filter::Filter(const std::string& in_fname) {
		//std::cerr << "Loading Filter from " << in_fname << '\n';

		//gzFile fp = gzopen(in_fname.c_str(),"rb");

		//std::ifstream infh(in_fname, std::ios::out | std::ios::binary);
		//gzread(fp, reinterpret_cast<char*>(&filter_size), sizeof(uint64_t));
		//gzread(fp, reinterpret_cast<char*>(&kmer_size), sizeof(uint64_t));
		//gzread(fp, reinterpret_cast<char*>(&hash_n), sizeof(uint64_t));

		//// Read the remaining data into a vector<uint8_t>
		//bytevec = std::vector<uint8_t>(filter_size);
		//gzread(fp, reinterpret_cast<char*>(bytevec.data()), filter_size);

		//std::cerr << "filter_size\t" << filter_size << '\n';
		//std::cerr << "kmer_size\t" << kmer_size << '\n';
		//std::cerr << "hash_n\t" << hash_n << '\n';

		//// Close the file
		//gzclose(fp);
	//}

	Filter::Filter(const std::string& in_fname) {
		std::cerr << "Loading Filter from " << in_fname << '\n';

		std::ifstream infh(in_fname, std::ios::out | std::ios::binary);
		infh.read(reinterpret_cast<char*>(&filter_size), sizeof(uint64_t));
		infh.read(reinterpret_cast<char*>(&kmer_size), sizeof(uint64_t));
		infh.read(reinterpret_cast<char*>(&hash_n), sizeof(uint64_t));

		// Read the remaining data into a vector<uint8_t>
		bytevec = std::vector<uint8_t>(filter_size);
		infh.read(reinterpret_cast<char*>(bytevec.data()), filter_size);

		std::cerr << "filter_size\t" << filter_size << '\n';
		std::cerr << "kmer_size\t" << kmer_size << '\n';
		std::cerr << "hash_n\t" << hash_n << '\n';

		// Close the file
		infh.close();
	}
	size_t Filter::setBitsCount() const {
		size_t count = 0;
		for (size_t i=0; i < size(); i++) {
			count += BitLookup::BIT_COUNT[bytevec[i]];
		}
		return count;
	}
}
