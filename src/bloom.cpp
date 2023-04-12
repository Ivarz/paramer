#include "bloom.h"
#include "bit_lookup.h"
#include <algorithm>
#include <cstdio>
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


	void Filter::writeRaw(const std::string& out_fname) const {
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

	void Filter::write(const std::string& out_fname) const {
		Gz::Writer gzwriter(out_fname);
		//gzFile fp = gzopen(out_fname.c_str(),"wb");

		uint64_t fsize = static_cast<uint64_t>(bytevec.size());
		uint64_t k = static_cast<uint64_t>(kmer_size);
		uint64_t h = static_cast<uint64_t>(hash_n);
		gzwriter.write(&fsize, sizeof(fsize));
		gzwriter.write(&k, sizeof(k));
		gzwriter.write(&h, sizeof(h));

		gzwriter.bufferedWrite(bytevec);
		//int max_bytes = std::numeric_limits<int>::max();
		//size_t offset = 0;

		//uint64_t written_bytes = 0;

		//while (written_bytes < filter_size) {
			//int buffer_size = std::min(static_cast<uint64_t>(max_bytes), filter_size - written_bytes);
			//int gzwrite_output = gzwrite(fp, (char*) &bytevec[offset], buffer_size*sizeof(bytevec.at(offset)));
			////std::cerr << buffer_size << '\t' << loaded_bytes << '\t' << filter_size << '\n';
			//if (gzwrite_output < 0) {
				////std::cerr << __FUNCTION__ << " Error " << gzerror(fp, &gzwrite_output)  << "\n";
				//gzclose(fp);
				////std::remove(out_fname);
				//return;
			//} else {
				//written_bytes += buffer_size;
				//offset = written_bytes;
			//}
		//}
		//std::cout << "written_bytes " << written_bytes << '\n';

		////outfh.write((char*) &out_fname[0], out_fname.size()*sizeof(char));
		////gzwrite(fp, (char*) &bytevec[0], filter_size*sizeof(bytevec.at(0)));

		//gzclose(fp);
	}

	std::optional<Filter> Filter::loadRaw(const std::string& in_fname) {
		uint64_t filter_size = 0;
		uint64_t kmer_size = 0;
		uint64_t hash_n = 0;

		std::ifstream infh(in_fname, std::ios::out | std::ios::binary);
		infh.read(reinterpret_cast<char*>(&filter_size), sizeof(uint64_t));
		infh.read(reinterpret_cast<char*>(&kmer_size), sizeof(uint64_t));
		infh.read(reinterpret_cast<char*>(&hash_n), sizeof(uint64_t));

		// Read the remaining data into a vector<uint8_t>
		Filter result = Filter(filter_size, kmer_size, hash_n);
		infh.read(reinterpret_cast<char*>(result.bytevec.data()), filter_size);

		// Close the file
		infh.close();
		return result;
	}

	std::optional<Filter> Filter::load(const std::string& in_fname) {

		uint64_t filter_size = 0;
		uint64_t kmer_size = 0;
		uint64_t hash_n = 0;

		Gz::Reader gz_reader(in_fname);

		gz_reader.read(&filter_size, sizeof(uint64_t));
		gz_reader.read(&kmer_size, sizeof(uint64_t));
		gz_reader.read(&hash_n, sizeof(uint64_t));

		// Read the remaining data into a vector<uint8_t>
		auto bufload = gz_reader.bufferedLoad(filter_size);
		if (bufload) {
			Filter result = Filter(0, 0, 0);
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
}
