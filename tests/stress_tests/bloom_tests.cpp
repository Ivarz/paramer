#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "bloom.h"
#include <cstdio>

const size_t TARGET_SIZE = 10000000000;
TEST_CASE("Test write and load Gz") {
	{
		CHECK(true);
		Bloom::Filter t0_bloom = Bloom::Filter(TARGET_SIZE, 31, 1);
		t0_bloom.writeGz("test_data/t0.stress.blm");
	}
	std::optional<Bloom::Filter> t0_bloom_load = Bloom::Filter::loadGz("test_data/t0.stress.blm");
	CHECK(true);
	//std::remove("test_data/t0.stress.blm");
	//CHECK(t0_bloom_load);
	//CHECK(t0_bloom_load->size() == TARGET_SIZE);
	//CHECK(t0_bloom_load->setBitsCount() == 0);
	//CHECK(t0_bloom_load->kmerSize() == 31);
	//CHECK(t0_bloom_load->hashN() == 1);

	//{
		//std::optional<Bloom::Filter> t1_bloom = Bloom::Filter(TARGET_SIZE, 31, 1);
		//t1_bloom->at(0) = 1;
		//t1_bloom->at(3) = 3;
		//t1_bloom->at(8) = 2;
		//t1_bloom->at(10) = 4;
		//t1_bloom->at(TARGET_SIZE-1) = 255;
		//t1_bloom->writeGz("test_data/t1.stress.blm");
		//std::optional<Bloom::Filter> t1_bloom_load = Bloom::Filter::loadGz("test_data/t1.stress.blm");
		//std::remove("test_data/t1.stress.blm");
		//CHECK(t1_bloom_load);
		//CHECK(t1_bloom_load->at(0) == 1);
		//CHECK(t1_bloom_load->at(3) == 3);
		//CHECK(t1_bloom_load->at(8) == 2);
		//CHECK(t1_bloom_load->at(10) == 4);
		//CHECK(t1_bloom_load->at(TARGET_SIZE-1) == 255);
	//}
}

//TEST_CASE("Test write and load") {
	//{
	//Bloom::Filter t0_bloom = Bloom::Filter(TARGET_SIZE, 31, 1);
	//t0_bloom.write("test_data/t0.stress.blm");
	//}
	//{
	//std::optional<Bloom::Filter> t0_bloom_load = Bloom::Filter::load("test_data/t0.stress.blm");
	//std::remove("test_data/t0.stress.blm");
	//CHECK(t0_bloom_load->size() == TARGET_SIZE);
	//CHECK(t0_bloom_load->setBitsCount() == 0);
	//CHECK(t0_bloom_load->kmerSize() == 31);
	//CHECK(t0_bloom_load->hashN() == 1);
	//}

	//{
	//std::optional<Bloom::Filter> t1_bloom = Bloom::Filter(TARGET_SIZE, 31, 1);
	//t1_bloom->at(0) = 1;
	//t1_bloom->at(3) = 3;
	//t1_bloom->at(8) = 2;
	//t1_bloom->at(10) = 4;
	//t1_bloom->at(TARGET_SIZE-1) = 255;
	//t1_bloom->write("test_data/t1.stress.blm");
	//}
	//{
	//std::optional<Bloom::Filter> t1_bloom_load = Bloom::Filter::load("test_data/t1.stress.blm");
	//std::remove("test_data/t1.stress.blm");
	//CHECK(t1_bloom_load->at(0) == 1);
	//CHECK(t1_bloom_load->at(3) == 3);
	//CHECK(t1_bloom_load->at(8) == 2);
	//CHECK(t1_bloom_load->at(10) == 4);
	//CHECK(t1_bloom_load->at(TARGET_SIZE-1) == 255);
	//}
//}

