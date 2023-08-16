#include "doctest.h"
#include "cmd.h"
#include <cstdio>
#include <vector>
#include "utils.h"
#include "fastx.h"
#include <fstream>


//test_data/test2.fa
//>seq1
//CCGATGCATCGATCGCCGACGACCACCCCCCCCCCCCCCCCCATGCATGATTCGACTGAG
//GTTTCCGTTGTAGGAGTGAACCCACTTGGCTTTGCGCCATAATTCCAATGAAAAACCTAT
//GCACTTTGTTTAGGGTACCATCAGGAATCTGAACCCTCAGATAGTGGGGATCCCGGGTAT
//>seq2
//AGACCTTTATCTGCGGTCCAACTTAGGCATAAACCTGCATGCTACCTTGTCAGACCCACT
//CTGCACGAAGTAAATATGGGATGCGTCCGACCTGGCTCCTGGCGTTCCACGCCGCCACGT
//GTTCGTTAACTGTTGATTGGTGGCACATAAGTAATACCATGGTCCCTGAAATTCGGCTCA
//>seq3
//GTTACTTCGAGCGTAATGTCTCAAATGGCGTAGAACGGCAATGACTGTTTGACACTAGGT
//GGTGTTCAGTTCGGTAACGGAGAGTCTGTGCGGCATTCTTATTAATACATTTGAAACGCG
//CCCAACTGACGCTAGGCAAGTCAGTGCAGGCTCCCGTGTTAGGATAAGGGTAAACATACA
//
//test_data/test2.k1.out.txt
//U       seq1    0       180     0:146
//U       seq2    0       180     1:146
//U       seq3    0       180     0:3 1:10 0:10 1:3 0:50 1:70
TEST_CASE("Test Cmd::Mask::run") {
	{
		int argn = 8;
		char* args[argn];
		args[0] = const_cast<char*>("paramer");
		args[1] = const_cast<char*>("mask");
		args[2] = const_cast<char*>("-k");
		args[3] = const_cast<char*>("test_data/test2.k1.out.txt");
		args[4] = const_cast<char*>("-f");
		args[5] = const_cast<char*>("test_data/test2.fa");
		args[6] = const_cast<char*>("-o");
		args[7] = const_cast<char*>("test_data/test2.mask.fa");
		Cmd::Mask::run(argn, args);

		Gz::Reader gzr("test_data/test2.mask.fa");
		std::optional<Fasta::Rec> rec1 = Fasta::nextRecord(gzr);
		std::optional<Fasta::Rec> rec2 = Fasta::nextRecord(gzr);
		std::optional<Fasta::Rec> rec3 = Fasta::nextRecord(gzr);
		CHECK(rec1->seq == "CCGATGCATCGATCGCCGACGACCACCCCCCCCCCCCCCCCCATGCATGATTCGACTGAGGTTTCCGTTGTAGGAGTGAACCCACTTGGCTTTGCGCCATAATTCCAATGAAAAACCTATGCACTTTGTTTAGGGTACCATCAGGAATCTGAACCCTCAGATAGTGGGGATCCCGGGTAT");
		CHECK(rec2->seq == "agacctttatctgcggtccaacttaggcataaacctgcatgctaccttgtcagacccactctgcacgaagtaaatatgggatgcgtccgacctggctcctggcgttccacgccgccacgtgttcgttaactgttgattggtggcacataagtaataccatggtccctgaaattcggctca");
		CHECK(rec3->seq == "GTTacttcgagcgtaatgtctcaaatggcgtagaacggcaatgactgtttgacactaggtGGTGTTCAGTTCGGTAacggagagtctgtgcggcattcttattaatacatttgaaacgcgcccaactgacgctaggcaagtcagtgcaggctcccgtgttaggataagggtaaacataca");
		std::remove("test_data/test2.mask.fa");
	}
	{
		int argn = 8;
		char* args[argn];
		args[0] = const_cast<char*>("paramer");
		args[1] = const_cast<char*>("mask");
		args[2] = const_cast<char*>("-k");
		args[3] = const_cast<char*>("test_data/test2.k1.out.txt.gz");
		args[4] = const_cast<char*>("-f");
		args[5] = const_cast<char*>("test_data/test2.fa.gz");
		args[6] = const_cast<char*>("-o");
		args[7] = const_cast<char*>("test_data/test2.mask.fa.gz");
		Cmd::Mask::run(argn, args);

		Gz::Reader gzr("test_data/test2.mask.fa.gz");
		std::optional<Fasta::Rec> rec1 = Fasta::nextRecord(gzr);
		std::optional<Fasta::Rec> rec2 = Fasta::nextRecord(gzr);
		std::optional<Fasta::Rec> rec3 = Fasta::nextRecord(gzr);
		CHECK(rec1->seq == "CCGATGCATCGATCGCCGACGACCACCCCCCCCCCCCCCCCCATGCATGATTCGACTGAGGTTTCCGTTGTAGGAGTGAACCCACTTGGCTTTGCGCCATAATTCCAATGAAAAACCTATGCACTTTGTTTAGGGTACCATCAGGAATCTGAACCCTCAGATAGTGGGGATCCCGGGTAT");
		CHECK(rec2->seq == "agacctttatctgcggtccaacttaggcataaacctgcatgctaccttgtcagacccactctgcacgaagtaaatatgggatgcgtccgacctggctcctggcgttccacgccgccacgtgttcgttaactgttgattggtggcacataagtaataccatggtccctgaaattcggctca");
		CHECK(rec3->seq == "GTTacttcgagcgtaatgtctcaaatggcgtagaacggcaatgactgtttgacactaggtGGTGTTCAGTTCGGTAacggagagtctgtgcggcattcttattaatacatttgaaacgcgcccaactgacgctaggcaagtcagtgcaggctcccgtgttaggataagggtaaacataca");
		std::remove("test_data/test2.mask.fa.gz");
	}
}

TEST_CASE("Test Cmd::BloomBuild::run") {
}
TEST_CASE("Test Cmd::StatsFasta::run") {
	{
		int argn = 6;
		char* args[argn];
		args[0] = const_cast<char*>("paramer");
		args[1] = const_cast<char*>("stats-fasta");
		args[2] = const_cast<char*>("-i");
		args[3] = const_cast<char*>("test_data/stats_t1.fa");
		args[4] = const_cast<char*>("-o");
		args[5] = const_cast<char*>("test_data/stats_t1.out.txt");
		Cmd::StatsFasta::run(argn, args);
		std::ifstream ifh("test_data/stats_t1.out.txt");
		std::string f1, f2, f3, f4;
		ifh >> f1 >> f2 >> f3 >> f4;
		CHECK(f1 == "seq_id");
		CHECK(f2 == "size");
		CHECK(f3 == "softmasked");
		CHECK(f4 == "hardmasked");
		ifh >> f1 >> f2 >> f3 >> f4;
		CHECK(f1 == "seq0");
		CHECK(f2 == "0");
		CHECK(f3 == "0");
		CHECK(f4 == "0");
		ifh >> f1 >> f2 >> f3 >> f4;
		CHECK(f1 == "seq1");
		CHECK(f2 == "12");
		CHECK(f3 == "0");
		CHECK(f4 == "0");
		ifh >> f1 >> f2 >> f3 >> f4;
		CHECK(f1 == "seq2");
		CHECK(f2 == "12");
		CHECK(f3 == "3");
		CHECK(f4 == "0");
		ifh >> f1 >> f2 >> f3 >> f4;
		CHECK(f1 == "seq3");
		CHECK(f2 == "12");
		CHECK(f3 == "0");
		CHECK(f4 == "2");
		ifh >> f1 >> f2 >> f3 >> f4;
		CHECK(f1 == "seq4");
		CHECK(f2 == "18");
		CHECK(f3 == "4");
		CHECK(f4 == "2");
		ifh >> f1 >> f2 >> f3 >> f4;
		CHECK(f1 == "seq5");
		CHECK(f2 == "28");
		CHECK(f3 == "12");
		CHECK(f4 == "4");
		ifh >> f1 >> f2 >> f3 >> f4;
		CHECK(f1 == "seq6");
		CHECK(f2 == "28");
		CHECK(f3 == "12");
		CHECK(f4 == "4");
		ifh.close();
		std::remove("test_data/stats_t1.out.txt");
	}
}
