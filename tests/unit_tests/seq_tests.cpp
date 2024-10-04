#include <string>
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "seq.h"
#include <iostream>

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

TEST_CASE("Testing Dna::unmask") {
	std::string t1 = "";
	Dna::unmask(t1,0,0);
	CHECK(t1 == "");

	std::string t2 = "atgc";
	Dna::unmask(t2,1,1);
	CHECK(t2 == "atgc");

	std::string t3 = "atgc";
	Dna::unmask(t3,1,2);
	CHECK(t3 == "aTgc");

	std::string t4 = "atgc";
	Dna::unmask(t4,1,3);
	CHECK(t4 == "aTGc");
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

TEST_CASE("Testing Dna::canonicalKmer") {
	std::string t0 = "";
	std::string t0_rc = Dna::revcom(t0);
	CHECK(Dna::canonicalKmer(t0) == t0);
	std::string t1 = "ATGCATCGATCG";
	std::string t1_rc = Dna::revcom(t1);
	CHECK(Dna::canonicalKmer(t1) == t1);
	std::string t2 = "GCATCGATCGATCGATGG";
	std::string t2_rc = Dna::revcom(t2);
	CHECK(Dna::canonicalKmer(t2) == t2_rc);
}
TEST_CASE("Testing Dna::getHashes") {
	std::string t1 = "AGTGCGTCGTCGTCGTCAGAGTGAAAACGTG";
	std::string t1_rc = Dna::revcom(t1);
	auto t1_hashes = Dna::getHashes(t1, 31, 1);
	auto t1_hashes_rc = Dna::getHashes(t1_rc, 31, 1);
	//std::cerr << t1 << '\t' << t1_hashes[0] << '\n';
	//std::cerr << t1_rc << '\t' << t1_hashes_rc[0] << '\n';
	CHECK(t1_hashes[0] == t1_hashes_rc[0]);
}

TEST_CASE("Testing Dna::getMinimizerHashes") {
	//std::vector<uint64_t> hashes = {2129312363048739376, 9967842831092148701, 4304134488848462692, 5819929052895302100, 12220444832806783367, 10118878782716752413, 6763474896827383532, 12220444832806783367, 10118878782716752413, 6763474896827383532, 12220444832806783367, 10118878782716752413, 6763474896827383532, 2737971792425662502, 11230458238545705290, 10883788240495922400, 2634314553780154417, 5922399782513233315, 1525707548702619386, 1194839728021048584, 5851829680807274824, 2555383599579742470, 1842717004517856612, 10829959314427631260, 11374462795461361796, 9484111255990750150, 5917157688684197953};
   
	{
		std::vector<uint64_t> true_result = 
		{
			2129312363048739376UL              
			, 4304134488848462692UL              
			, 4304134488848462692UL              
			, 5819929052895302100UL              
			, 6763474896827383532UL              
			, 6763474896827383532UL              
			, 6763474896827383532UL              
			, 6763474896827383532UL              
			, 6763474896827383532UL              
			, 6763474896827383532UL              
			, 2737971792425662502UL              
			, 2737971792425662502UL              
			, 2737971792425662502UL              
			, 2634314553780154417UL              
			, 2634314553780154417UL              
			, 1525707548702619386UL              
			, 1194839728021048584UL              
			, 1194839728021048584UL              
			, 1194839728021048584UL              
			, 1194839728021048584UL              
			, 1842717004517856612UL              
			, 1842717004517856612UL              
			, 1842717004517856612UL              
			, 5917157688684197953UL              
		};
		std::vector<uint64_t> minimizers = Dna::getMinimizerHashes("AGTGCGTCGTCGTCGTCAGAGTGAAAACGTG", 5, 1, 8);
		for (size_t i=0; i < true_result.size(); i++) {
			CHECK( minimizers[i] == true_result[i]);
		}
	}
	{
		std::vector<uint64_t> true_result = 
		{
			2129312363048739376UL
			, 2183400668632071988UL
			, 7717741254769009360UL
			, 4304134488848462692UL
			, 2183400668632071988UL
			, 130905002315143641UL
			, 4304134488848462692UL
			, 2183400668632071988UL
			, 130905002315143641UL
			, 5819929052895302100UL
			, 7525786875638661037UL
			, 130905002315143641UL
			, 6763474896827383532UL
			, 7525786875638661037UL
			, 130905002315143641UL
			, 6763474896827383532UL
			, 7525786875638661037UL
			, 130905002315143641UL
			, 6763474896827383532UL
			, 7525786875638661037UL
			, 130905002315143641UL
			, 6763474896827383532UL
			, 7525786875638661037UL
			, 130905002315143641UL
			, 6763474896827383532UL
			, 7525786875638661037UL
			, 130905002315143641UL
			, 6763474896827383532UL
			, 7525786875638661037UL
			, 130905002315143641UL
			, 2737971792425662502UL
			, 592103693697923601UL
			, 130905002315143641UL
			, 2737971792425662502UL
			, 592103693697923601UL
			, 1384494391531451430UL
			, 2737971792425662502UL
			, 592103693697923601UL
			, 5682106283415464531UL
			, 2634314553780154417UL
			, 343377397220436461UL
			, 6665969643078428448UL
			, 2634314553780154417UL
			, 343377397220436461UL
			, 6665969643078428448UL
			, 1525707548702619386UL
			, 343377397220436461UL
			, 4828501057582604783UL
			, 1194839728021048584UL
			, 343377397220436461UL
			, 4828501057582604783UL
			, 1194839728021048584UL
			, 5271014183005424866UL
			, 4828501057582604783UL
			, 1194839728021048584UL
			, 5271014183005424866UL
			, 4828501057582604783UL
			, 1194839728021048584UL
			, 5271014183005424866UL
			, 4845804533368155100UL
			, 1842717004517856612UL
			, 5271014183005424866UL
			, 4845804533368155100UL
			, 1842717004517856612UL
			, 5482292490891933263UL
			, 4845804533368155100UL
			, 1842717004517856612UL
			, 5482292490891933263UL
			, 9885902586050233954UL
			, 5917157688684197953UL
			, 5482292490891933263UL
			, 9587070619738185232UL

		};
		std::vector<uint64_t> minimizers = Dna::getMinimizerHashes("AGTGCGTCGTCGTCGTCAGAGTGAAAACGTG", 5, 3, 8);
		for (size_t i=0; i < true_result.size(); i++) {
			CHECK( minimizers[i] == true_result[i]);
		}
	}
}

TEST_CASE("Testing Dna::shannon") {
	CHECK(std::abs(Dna::shannon("") - 0.0) < 0.000001 );
	CHECK(std::abs(Dna::shannon("ATGATGATGATGATGATGATGATG") - 1.0930808359255935) < 0.000001 );
	CHECK(std::abs(Dna::shannon("GTCTGTCTGGTAGCTTACGATTTCCCGTGGT") - 3.0599691653497123) < 0.000001 );
}

TEST_CASE("Test Dna::maskingStats") {
	{
		std::string test_case = "";
		Dna::MaskingStats res = Dna::maskingStats(test_case);
		CHECK(res.size == 0);
		CHECK(res.softmasked == 0);
		CHECK(res.hardmasked == 0);
	}
	{
		std::string test_case = "A";
		Dna::MaskingStats res = Dna::maskingStats(test_case);
		CHECK(res.size == 1);
		CHECK(res.softmasked == 0);
		CHECK(res.hardmasked == 0);
	}
	{
		std::string test_case = "a";
		Dna::MaskingStats res = Dna::maskingStats(test_case);
		CHECK(res.size == 1);
		CHECK(res.softmasked == 1);
		CHECK(res.hardmasked == 0);
	}
	{
		std::string test_case = "N";
		Dna::MaskingStats res = Dna::maskingStats(test_case);
		CHECK(res.size == 1);
		CHECK(res.softmasked == 0);
		CHECK(res.hardmasked == 1);
	}
	{
		std::string test_case = "n";
		Dna::MaskingStats res = Dna::maskingStats(test_case);
		CHECK(res.size == 1);
		CHECK(res.softmasked == 0);
		CHECK(res.hardmasked == 1);
	}
	{
		std::string test_case = "atgcn";
		Dna::MaskingStats res = Dna::maskingStats(test_case);
		CHECK(res.size == 5);
		CHECK(res.softmasked == 4);
		CHECK(res.hardmasked == 1);
	}
	{
		std::string test_case = "ATGCatgcNn";
		Dna::MaskingStats res = Dna::maskingStats(test_case);
		CHECK(res.size == 10);
		CHECK(res.softmasked == 4);
		CHECK(res.hardmasked == 2);
	}
}
