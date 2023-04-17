#include "doctest.h"
#include "bloom.h"
#include "bit_lookup.h"
#include <cstdio>

TEST_CASE("Test constructor") {
	Bloom::Filter t1 = Bloom::Filter(1000, 31, 31, 3);
	CHECK(t1.size() == 1000);
	CHECK(t1.kmerSize() == 31);
	CHECK(t1.hashN() == 3);
}

TEST_CASE("Test bit lookup table") {
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

TEST_CASE("Test setBitsCount") {
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

TEST_CASE("Test addSeq") {
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

TEST_CASE("Test write and load") {
	Bloom::Filter t0_bloom = Bloom::Filter(1000, 31, 31, 1);
	t0_bloom.write("test_data/t0.blm");
	std::optional<Bloom::Filter> t0_bloom_load = Bloom::Filter::load("test_data/t0.blm");
	std::remove("test_data/t0.blm");
	CHECK(t0_bloom_load->size() == 1000);
	CHECK(t0_bloom_load->setBitsCount() == 0);
	CHECK(t0_bloom_load->kmerSize() == 31);
	CHECK(t0_bloom_load->hashN() == 1);

	Bloom::Filter t1_bloom = Bloom::Filter(1000, 31, 31, 1);
	t1_bloom.at(0) = 1;
	t1_bloom.at(3) = 3;
	t1_bloom.at(8) = 2;
	t1_bloom.at(10) = 4;
	t1_bloom.at(999) = 255;
	t1_bloom.write("test_data/t1.blm");
	std::optional<Bloom::Filter> t1_bloom_load = Bloom::Filter::load("test_data/t1.blm");
	std::remove("test_data/t1.blm");
	CHECK(t1_bloom_load->at(0) == 1);
	CHECK(t1_bloom_load->at(3) == 3);
	CHECK(t1_bloom_load->at(8) == 2);
	CHECK(t1_bloom_load->at(10) == 4);
	CHECK(t1_bloom_load->at(999) == 255);
}

TEST_CASE("Test write and load Raw") {
	Bloom::Filter t0_bloom = Bloom::Filter(1000, 31, 31, 1);
	t0_bloom.writeRaw("test_data/t0.blm");
	std::optional<Bloom::Filter> t0_bloom_load = Bloom::Filter::loadRaw("test_data/t0.blm");
	std::remove("test_data/t0.blm");
	CHECK(t0_bloom_load->size() == 1000);
	CHECK(t0_bloom_load->setBitsCount() == 0);
	CHECK(t0_bloom_load->kmerSize() == 31);
	CHECK(t0_bloom_load->hashN() == 1);

	Bloom::Filter t1_bloom = Bloom::Filter(1000, 31, 31, 1);
	t1_bloom.at(0) = 1;
	t1_bloom.at(3) = 3;
	t1_bloom.at(8) = 2;
	t1_bloom.at(10) = 4;
	t1_bloom.at(999) = 255;
	t1_bloom.writeRaw("test_data/t1.blm");
	std::optional<Bloom::Filter> t1_bloom_load = Bloom::Filter::loadRaw("test_data/t1.blm");
	std::remove("test_data/t1.blm");
	CHECK(t1_bloom_load->at(0) == 1);
	CHECK(t1_bloom_load->at(3) == 3);
	CHECK(t1_bloom_load->at(8) == 2);
	CHECK(t1_bloom_load->at(10) == 4);
	CHECK(t1_bloom_load->at(999) == 255);
}

TEST_CASE ("Test extendSeq") {
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
		std::string seq2 = "CACCTGTATTTTTCGGTTCTGCACTGACAAATTTTGGTGTGGAAACATTTTTAAAACATCGTTACGCACAGATGACATATTCTCCATTGCCACAATTTTTACAGATGACATATTCTCCATTGCCA";
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
