#include "doctest.h"
#include "bloom.h"

TEST_CASE("Test constructor") {
	Bloom::Filter t1 = Bloom::Filter(1000, 31, 3);
	CHECK(t1.size() == 1000);
	CHECK(t1.kmerSize() == 31);
	CHECK(t1.hashN() == 3);
}

TEST_CASE("Test setBitsCount") {
	std::string t1_seq = "AGTGCGTCGTCGTCGTCAGAGTGAAAACGTG";
	Bloom::Filter t1_bloom = Bloom::Filter(1000, 31, 1);
	t1_bloom.addSeq(t1_seq);
	CHECK(t1_bloom.setBitsCount() == 1);

	std::string t2_seq = "AGTGCGTCGTCGTCGTCAGAGTGAAAACGTG";
	Bloom::Filter t2_bloom = Bloom::Filter(1000, 31, 3);
	t2_bloom.addSeq(t2_seq);
	CHECK(t2_bloom.setBitsCount() == 3);
}

TEST_CASE("Test addSeq") {
	std::string t0_seq = "";
	Bloom::Filter t0_bloom = Bloom::Filter(1000, 31, 1);
	t0_bloom.addSeq(t0_seq);
	CHECK(t0_bloom.setBitsCount() == 0);

	std::string t1_seq = "AGTGCGTCGTCGTCGTCAGAGTGAAAACGTG";
	uint64_t t1_hash = 6248067913986390878;
	Bloom::Filter t1_bloom = Bloom::Filter(1000, 31, 1);
	uint64_t t1_bit_idx = t1_hash % (t1_bloom.size()*8); // 7024th bit
	uint64_t t1_byte_idx = t1_bit_idx / 8; //859th byte
	uint8_t byte_value = 1 << (t1_bit_idx % 8); // 859th bytes value should be 64
	t1_bloom.addSeq(t1_seq);
	CHECK(t1_bloom.at(t1_byte_idx) == byte_value);
}
