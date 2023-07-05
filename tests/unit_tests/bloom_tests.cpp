#include "doctest.h"
#include "bloom.h"
#include "bit_lookup.h"
#include <cstdio>

TEST_CASE("Test Bloom::Filter constructor") {
	Bloom::Filter t1 = Bloom::Filter(1000, 31, 31, 3);
	CHECK(t1.size() == 1000);
	CHECK(t1.kmerSize() == 31);
	CHECK(t1.hashN() == 3);
}

TEST_CASE("Test BitLookup::BIT_COUNT table") {
	for (size_t i = 1; i < 256; i++){
		uint8_t table_value = BitLookup::BIT_COUNT[i];
		unsigned int calculated_value = 0;
		size_t n = i;
		while (n > 0) {
			n = n & (n-1);
			calculated_value++;
		}
		CHECK(table_value == calculated_value);
	}
}

TEST_CASE("Test Bloom::Filter::setBitsCount") {
	Bloom::Filter t0_bloom = Bloom::Filter(1000, 31, 31, 1);
	CHECK(t0_bloom.setBitsCount() == 0);

	std::string t1_seq = "AGTGCGTCGTCGTCGTCAGAGTGAAAACGTG";
	Bloom::Filter t1_bloom = Bloom::Filter(1000, 31, 31, 1);
	t1_bloom.addSeq(t1_seq);
	CHECK(t1_bloom.setBitsCount() == 1);

	std::string t2_seq = "AGTGCGTCGTCGTCGTCAGAGTGAAAACGTG";
	Bloom::Filter t2_bloom = Bloom::Filter(1000, 31, 31, 3);
	t2_bloom.addSeq(t2_seq);
	CHECK(t2_bloom.setBitsCount() == 3);
}

TEST_CASE("Test Bloom::Filter::addSeq") {
	std::string t0_seq = "";
	Bloom::Filter t0_bloom = Bloom::Filter(1000, 31, 31, 1);
	t0_bloom.addSeq(t0_seq);
	CHECK(t0_bloom.setBitsCount() == 0);

	std::string t1_seq = "AGTGCGTCGTCGTCGTCAGAGTGAAAACGTG";
	uint64_t t1_hash = 6248067913986390878;
	Bloom::Filter t1_bloom = Bloom::Filter(1000, 31, 31, 1);
	uint64_t t1_bit_idx = t1_hash % (t1_bloom.size()*8); // 7024th bit
	uint64_t t1_byte_idx = t1_bit_idx / 8; //859th byte
	uint8_t byte_value = 1 << (t1_bit_idx % 8); // 859th bytes value should be 64
	t1_bloom.addSeq(t1_seq);
	CHECK(t1_bloom.at(t1_byte_idx) == byte_value);
}

TEST_CASE("Test Bloom::Filter::searchSeq") {
	{
		std::string seq = "AGTGCGTCGTCGTCGTCAGAGTGAAAACGTGCGCATGACTGACTGACTGACGTACAGGAA";
		{
			Bloom::Filter bloom = Bloom::Filter(1000, 31, 35, 3);
			bloom.addSeq(seq);
			size_t result = bloom.searchSeq(seq);
			CHECK(result == 30);
			bloom.writeRaw("test_data/t1.blm");
		}
		{
			std::optional<Bloom::Filter> bloom_load = Bloom::Filter::loadPointer("test_data/t1.blm");
			size_t result = bloom_load->searchSeq(seq);
			CHECK(result == 30);
			std::remove("test_data/t1.blm");
		}
	}
}

TEST_CASE("Test Bloom::Filter::seekSeq") {
	{
		std::string seq = "AGTGCGTCGTCGTCGTCAGAGTGAAAACGTGCGCATGACTGACTGACTGACGTACAGGAA";
		{
			Bloom::Filter bloom = Bloom::Filter(1000, 31, 31, 3);
			bloom.addSeq(seq);
			bloom.writeRaw("test_data/t1.blm");
		}
		{
			std::optional<Bloom::Filter> bloom_load = Bloom::Filter::loadPointer("test_data/t1.blm");
			size_t result = bloom_load->seekSeq(seq);
			CHECK(result == 30);
			std::remove("test_data/t1.blm");
		}
	}
}

TEST_CASE("Test Bloom::Filter::searchMinimizers") {
	{
		std::string seq = "AGTGCGTCGTCGTCGTCAGAGTGAAAACGTGCGCATGACTGACTGACTGACGTACAGGAA";
		{
			Bloom::Filter bloom = Bloom::Filter(1000, 31, 35, 3);
			bloom.addSeq(seq);
			size_t result = bloom.searchMinimizers(seq);
			CHECK(result == 26);
			bloom.writeRaw("test_data/t1.blm");
		}
		{
			std::optional<Bloom::Filter> bloom_load = Bloom::Filter::loadPointer("test_data/t1.blm");
			size_t result = bloom_load->searchMinimizers(seq);
			CHECK(result == 26);
			std::remove("test_data/t1.blm");
		}
	}
}


TEST_CASE("Test Bloom::Filter::searchFastqPair") {
	{
		Gz::Reader reader1 = Gz::Reader("test_data/test.sub.1.fq.gz");
		Gz::Reader reader2 = Gz::Reader("test_data/test.sub.2.fq.gz");
		std::optional<Fastq::Pair> rec_pair = Fastq::nextRecordPair(reader1, reader2);
		std::string seq = "AATTTGACATGGATCTTGTATCAAAGGGAGAACTTTCACCTGTATTTTTCGGTTCTGCAC"; //subsequence from test.sub.1.fq.gz
		Bloom::Filter bloom = Bloom::Filter(1000, 31, 31, 3);
		bloom.addSeq(seq);
		size_t result = bloom.searchFastqPair(*rec_pair);
		CHECK(result == 30);
	}
	{
		Gz::Reader reader1 = Gz::Reader("test_data/test.sub.1.fq.gz");
		Gz::Reader reader2 = Gz::Reader("test_data/test.sub.2.fq.gz");
		std::optional<Fastq::Pair> rec_pair = Fastq::nextRecordPair(reader1, reader2);
		std::string seq = "ATTTACCGGAACAGATACGCATAAATGCGATCCTGTCCCTGTGGTTTTTATTCATATTTG"; //subsequence from test.sub.2.fq.gz
		Bloom::Filter bloom = Bloom::Filter(1000, 31, 31, 3);
		bloom.addSeq(seq);
		size_t result = bloom.searchFastqPair(*rec_pair);
		CHECK(result == 30);
	}
}

TEST_CASE("Test Bloom::Filter::writeGz and Bloom::Filter::loadGz") {
	Bloom::Filter t0_bloom = Bloom::Filter(1000, 31, 31, 1);
	t0_bloom.writeGz("test_data/t0.blm");
	std::optional<Bloom::Filter> t0_bloom_load = Bloom::Filter::loadGz("test_data/t0.blm");
	std::remove("test_data/t0.blm");
	CHECK(t0_bloom_load->size() == 1000);
	CHECK(t0_bloom_load->setBitsCount() == 0);
	CHECK(t0_bloom_load->kmerSize() == 31);
	CHECK(t0_bloom_load->hashN() == 1);

	Bloom::Filter t1_bloom = Bloom::Filter(1000, 31, 31, 1);
	t1_bloom.atRef(0) = 1;
	t1_bloom.atRef(3) = 3;
	t1_bloom.atRef(8) = 2;
	t1_bloom.atRef(10) = 4;
	t1_bloom.atRef(999) = 255;
	t1_bloom.writeGz("test_data/t1.blm");
	std::optional<Bloom::Filter> t1_bloom_load = Bloom::Filter::loadGz("test_data/t1.blm");
	std::remove("test_data/t1.blm");
	CHECK(t1_bloom_load->at(0) == 1);
	CHECK(t1_bloom_load->at(3) == 3);
	CHECK(t1_bloom_load->at(8) == 2);
	CHECK(t1_bloom_load->at(10) == 4);
	CHECK(t1_bloom_load->at(999) == 255);
}

TEST_CASE("Test Bloom::Filter::writeRaw and Bloom::Filter::loadRaw") {
	Bloom::Filter t0_bloom = Bloom::Filter(1000, 31, 31, 1);
	t0_bloom.writeRaw("test_data/t0.blm");
	std::optional<Bloom::Filter> t0_bloom_load = Bloom::Filter::loadRaw("test_data/t0.blm");
	std::remove("test_data/t0.blm");
	CHECK(t0_bloom_load->size() == 1000);
	CHECK(t0_bloom_load->setBitsCount() == 0);
	CHECK(t0_bloom_load->kmerSize() == 31);
	CHECK(t0_bloom_load->hashN() == 1);

	Bloom::Filter t1_bloom = Bloom::Filter(1000, 31, 31, 1);
	t1_bloom.atRef(0) = 1;
	t1_bloom.atRef(3) = 3;
	t1_bloom.atRef(8) = 2;
	t1_bloom.atRef(10) = 4;
	t1_bloom.atRef(999) = 255;
	t1_bloom.writeRaw("test_data/t1.blm");
	std::optional<Bloom::Filter> t1_bloom_load = Bloom::Filter::loadRaw("test_data/t1.blm");
	std::remove("test_data/t1.blm");
	CHECK(t1_bloom_load->at(0) == 1);
	CHECK(t1_bloom_load->at(3) == 3);
	CHECK(t1_bloom_load->at(8) == 2);
	CHECK(t1_bloom_load->at(10) == 4);
	CHECK(t1_bloom_load->at(999) == 255);
}

TEST_CASE("Test Bloom::Filter::write and Bloom::Filter::load") {
	{
		Bloom::Compression cmpr = Bloom::Compression::RAW;
		Bloom::Filter t0_bloom = Bloom::Filter(1000, 31, 31, 1);
		t0_bloom.write("test_data/t0.blm", cmpr);
		std::optional<Bloom::Filter> t0_bloom_load = Bloom::Filter::load("test_data/t0.blm");
		std::remove("test_data/t0.blm");
		CHECK(t0_bloom_load->size() == 1000);
		CHECK(t0_bloom_load->setBitsCount() == 0);
		CHECK(t0_bloom_load->kmerSize() == 31);
		CHECK(t0_bloom_load->hashN() == 1);

		Bloom::Filter t1_bloom = Bloom::Filter(1000, 31, 31, 1);
		t1_bloom.atRef(0) = 1;
		t1_bloom.atRef(3) = 3;
		t1_bloom.atRef(8) = 2;
		t1_bloom.atRef(10) = 4;
		t1_bloom.atRef(999) = 255;
		t1_bloom.write("test_data/t1.blm", cmpr);
		std::optional<Bloom::Filter> t1_bloom_load = Bloom::Filter::load("test_data/t1.blm");
		std::remove("test_data/t1.blm");
		CHECK(t1_bloom_load->at(0) == 1);
		CHECK(t1_bloom_load->at(3) == 3);
		CHECK(t1_bloom_load->at(8) == 2);
		CHECK(t1_bloom_load->at(10) == 4);
		CHECK(t1_bloom_load->at(999) == 255);
	}
	{
		Bloom::Compression cmpr = Bloom::Compression::GZ;
		Bloom::Filter t0_bloom = Bloom::Filter(1000, 31, 31, 1);
		t0_bloom.write("test_data/t0.blm", cmpr);
		std::optional<Bloom::Filter> t0_bloom_load = Bloom::Filter::load("test_data/t0.blm");
		std::remove("test_data/t0.blm");
		CHECK(t0_bloom_load->size() == 1000);
		CHECK(t0_bloom_load->setBitsCount() == 0);
		CHECK(t0_bloom_load->kmerSize() == 31);
		CHECK(t0_bloom_load->hashN() == 1);

		Bloom::Filter t1_bloom = Bloom::Filter(1000, 31, 31, 1);
		t1_bloom.atRef(0) = 1;
		t1_bloom.atRef(3) = 3;
		t1_bloom.atRef(8) = 2;
		t1_bloom.atRef(10) = 4;
		t1_bloom.atRef(999) = 255;
		t1_bloom.write("test_data/t1.blm", cmpr);
		std::optional<Bloom::Filter> t1_bloom_load = Bloom::Filter::load("test_data/t1.blm");
		std::remove("test_data/t1.blm");
		CHECK(t1_bloom_load->at(0) == 1);
		CHECK(t1_bloom_load->at(3) == 3);
		CHECK(t1_bloom_load->at(8) == 2);
		CHECK(t1_bloom_load->at(10) == 4);
		CHECK(t1_bloom_load->at(999) == 255);
	}
}

TEST_CASE ("Test Bloom::Filter::extendSeq") {
	{
		Bloom::Filter t1_bloom = Bloom::Filter(1000, 31, 31, 3);
		std::string t1_seq1 = "GAACTCTTAGACGGTGCAAGCGCAGAATTTGACATGGATCTTGTATCAAAGGGAGAACTTTCACCTGTATTTTTCGGTTCTGCACTGACAAATTTTGGTGTGGAAACATTTTTAAAACATTTTTTACAGATGACATATTCTCCATTGCCA";
		std::string t1_seq2 = "GAACTCTTAGACGGTGCAAGCGCAGAATTTGACATGGATCTTGTATCAAAGGGAGAACTTTCACCTGTATTTTTCGGTTCTGCACTGACAAATTTTGGTGTGGAAACATTTTTAAAGCGTTACGCGTACGATCGATCAGTCATGACTG";
		t1_bloom.addSeq(t1_seq1);
		t1_bloom.addSeq(t1_seq2);
		auto result = t1_bloom.extendSeq("ATCTTGTATCAAAGGGAGAACTTTCACCTGTATTTTTCGGTTCTGCACTGACAAATTTTGGTGTGGAAACATTTTTAAA");
		std::sort(result.begin(), result.end());
		CHECK(result.size() == 2);
		CHECK(result.at(0) == t1_seq1);
		CHECK(result.at(1) == t1_seq2);
	}

	{
		Bloom::Filter bloom = Bloom::Filter(1000, 31, 31, 3);
		std::string seq1 = "CACCTGTATTTTTCGGTTCTGCACTGACAAATTTTGGTGTGGAAACATTTTTAAAACAATTTTTACAGATGACATATTCTCCATTGCCA";
		std::string seq2 = "CACCTGTATTTTTCGGTTCTGCACTGACAAATTTTGGTGTGGAAACATTTTTAAAACATCGTCAATCGTACGATCGACCTCTCTCATGCTAACTGACTGACTGACTGCATG";
		bloom.addSeq(seq1);
		bloom.addSeq(seq2);
		auto result = bloom.extendSeq("CACCTGTATTTTTCGGTTCTGCACTGACAAATTTTGGTGTGGAAACATTTTTAAA");
		std::sort(result.begin(), result.end());
		CHECK(result.size() == 2);
		CHECK(result.at(0) == seq1);
		CHECK(result.at(1) == seq2);
	}
	//std::string t1_seq3 = "CGTACGTCGTCGCTGTGTCAGCAGAGAGGGGCCCGGCCGAGTAGAATCTTGTATCAAAGGGAGAACTTTCACCTGTATTTTTCGGTTCTGCACTGACAAATTTTGGTGTGGAAACATTTTTAAAACATTTTTTACAGATGACATATTCTCCATTGCCA";
	//t1_bloom.addSeq(t1_seq3);
	//auto result2 = t1_bloom.extendSeq("ATCTTGTATCAAAGGGAGAACTTTCACCTGTATTTTTCGGTTCTGCACTGACAAATTTTGGTGTGGAAACATTTTTAAA");
	//auto result = t1_bloom.extendSeq("ATCTTGTATCAAAGGGAGAACTTTCACCTGTATTTTTCGGTTCTGCACTGACAAATTTTGGTGTGGAAACATTTTTAAA");

}
