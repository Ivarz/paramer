#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "seq.h"

TEST_CASE("Testing Dna::revcom") {
	CHECK(Dna::revcom("") == "");
	CHECK(Dna::revcom("A") == "T");
	CHECK(Dna::revcom("T") == "A");
	CHECK(Dna::revcom("G") == "C");
	CHECK(Dna::revcom("C") == "G");
	CHECK(Dna::revcom("a") == "t");
	CHECK(Dna::revcom("t") == "a");
	CHECK(Dna::revcom("g") == "c");
	CHECK(Dna::revcom("c") == "g");
	CHECK(Dna::revcom("N") == "N");
	CHECK(Dna::revcom("n") == "n");
	CHECK(Dna::revcom("ATGC") == "GCAT");
	CHECK(Dna::revcom("atgc") == "gcat");
	CHECK(Dna::revcom("ATGCN") == "NGCAT");
	CHECK(Dna::revcom("atgcn") == "ngcat");
}

TEST_CASE("Testing Dna::softmask") {
	std::string t1 = "";
	Dna::softmask(t1,0,0);
	CHECK(t1 == "");

	std::string t2 = "ATGC";
	Dna::softmask(t2,1,1);
	CHECK(t2 == "ATGC");

	std::string t3 = "ATGC";
	Dna::softmask(t3,1,2);
	CHECK(t3 == "AtGC");

	std::string t4 = "ATGC";
	Dna::softmask(t4,1,3);
	CHECK(t4 == "AtgC");
}

TEST_CASE("Testing Dna::nextToggleMaskedRegion") {

	auto t1 = Dna::nextToggleMaskedRegion("", 0);
	CHECK(t1 == std::pair<size_t,size_t>(0,0));

	auto t2 = Dna::nextToggleMaskedRegion("", 123);
	CHECK(t2 == std::pair<size_t,size_t>(0,0));

	auto t3 = Dna::nextToggleMaskedRegion("A", 0);
	CHECK(t3 == std::pair<size_t,size_t>(0,1));

	auto t4 = Dna::nextToggleMaskedRegion("a", 0);
	CHECK(t4 == std::pair<size_t,size_t>(0,1));

	auto t5 = Dna::nextToggleMaskedRegion("Aa", 0);
	CHECK(t5 == std::pair<size_t,size_t>(0,1));

	auto t6 = Dna::nextToggleMaskedRegion("Aa", 1);
	CHECK(t6 == std::pair<size_t,size_t>(1,2));

	auto t7 = Dna::nextToggleMaskedRegion("aA", 0);
	CHECK(t7 == std::pair<size_t,size_t>(0,1));

	auto t8 = Dna::nextToggleMaskedRegion("aA", 1);
	CHECK(t8 == std::pair<size_t,size_t>(1,2));

	auto t9 = Dna::nextToggleMaskedRegion("ATGC", 0);
	CHECK(t9 == std::pair<size_t,size_t>(0,4));

	auto t10 = Dna::nextToggleMaskedRegion("atgc", 0);
	CHECK(t10 == std::pair<size_t,size_t>(0,4));

	auto t11 = Dna::nextToggleMaskedRegion("ATGCatgc", 0);
	CHECK(t11 == std::pair<size_t,size_t>(0,4));

	auto t12 = Dna::nextToggleMaskedRegion("ATGCatgc", 4);
	CHECK(t12 == std::pair<size_t,size_t>(4,8));

	auto t13 = Dna::nextToggleMaskedRegion("atgcATGC", 0);
	CHECK(t13 == std::pair<size_t,size_t>(0,4));

	auto t14 = Dna::nextToggleMaskedRegion("atgcATGC", 4);
	CHECK(t14 == std::pair<size_t,size_t>(4,8));
}

TEST_CASE("Testing Dna::splitOnMask") {
	std::vector<std::string> empty_res{};
	CHECK(Dna::splitOnMask("") == empty_res);
	CHECK(Dna::splitOnMask("atgc") == empty_res);

	std::vector<std::string> t2 = Dna::splitOnMask("ATGC");
	CHECK(t2[0] == "ATGC");
	CHECK(t2.size() == 1);

	std::vector<std::string> t3 = Dna::splitOnMask("ATaGC");
	CHECK(t3[0] == "AT");
	CHECK(t3[1] == "GC");
	CHECK(t3.size() == 2);

	std::vector<std::string> t4 = Dna::splitOnMask("ATtGC");
	CHECK(t4[0] == "AT");
	CHECK(t4[1] == "GC");
	CHECK(t4.size() == 2);

	std::vector<std::string> t5 = Dna::splitOnMask("ATgGC");
	CHECK(t5[0] == "AT");
	CHECK(t5[1] == "GC");
	CHECK(t5.size() == 2);

	std::vector<std::string> t6 = Dna::splitOnMask("ATcGC");
	CHECK(t6[0] == "AT");
	CHECK(t6[1] == "GC");
	CHECK(t6.size() == 2);

	std::vector<std::string> t7 = Dna::splitOnMask("ATNGC");
	CHECK(t7[0] == "AT");
	CHECK(t7[1] == "GC");
	CHECK(t7.size() == 2);

	std::vector<std::string> t8 = Dna::splitOnMask("ATnGC");
	CHECK(t8[0] == "AT");
	CHECK(t8[1] == "GC");
	CHECK(t8.size() == 2);

}

