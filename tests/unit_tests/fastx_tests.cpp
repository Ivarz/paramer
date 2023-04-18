#include "doctest.h"
#include "fastx.h"
#include "kraken2.h"

TEST_CASE("Testing Fasta::nextRecord") {
	{
		Gz::Reader reader("test_data/empty.fa");
		std::optional<Fasta::Rec> rec = Fasta::nextRecord(reader);
		CHECK(!rec);
	}
	{
		Gz::Reader reader("test_data/t1.fa.gz");
		std::optional<Fasta::Rec> rec1 = Fasta::nextRecord(reader);
		std::optional<Fasta::Rec> rec2 = Fasta::nextRecord(reader);
		std::optional<Fasta::Rec> rec3 = Fasta::nextRecord(reader);
		std::optional<Fasta::Rec> rec4 = Fasta::nextRecord(reader);
		std::optional<Fasta::Rec> rec_end = Fasta::nextRecord(reader);
		CHECK(rec1->seq_id == "EVEC_sCAffold0000001 lenGTh=176877");
		CHECK(rec1->seq.size() == 176877);
		CHECK(rec2->seq_id == "EVEC_sCAffold0000002 lenGTh=165254");
		CHECK(rec2->seq.size() == 165254);
		CHECK(rec3->seq_id == "EVEC_sCAffold0000003 lenGTh=155333");
		CHECK(rec3->seq.size() == 155333);
		CHECK(rec4->seq_id == "EVEC_sCAffold0000004 lenGTh=147707");
		CHECK(rec4->seq.size() == 147707);
		CHECK(!rec_end);
	}
}

TEST_CASE("Testing Fasta::softmaskWithKraken2") {
	{
		Gz::Reader k2_reader("test_data/t3.k2out.txt");
		Gz::Reader fa_reader("test_data/t1.fa.gz");
		std::optional<Kraken2::Rec> krec1 = Kraken2::nextRecord(k2_reader);
		std::optional<Fasta::Rec> farec1 = Fasta::nextRecord(fa_reader);
		CHECK(farec1->seq.substr(0,60) == "GCAAAAAAGTAAACAAAAAGCAAAACAAACAGAAAAAAAAGAAAAAAATATTTGTTTTTT");
		CHECK(farec1->seq.substr(1729,70) == "AAAGACTGACTGTTTAATTGTATTTTAACGTATCCCACTTCGAACCTCTAACAACTGCGGGGCGATATCG");
		farec1->softmaskWithKraken2(*krec1, 35);
		CHECK(farec1->seq.substr(0,60) == "gcaaaaaagtaaacaaaaagcaaaacaaacagaaaaaaaagaaaaaaataTTTGTTTTTT");
		CHECK(farec1->seq.substr(1729,70) == "AAAGACTGACtgtttaattgtattttaacgtatcccacttcgaacctctaacaactgcggGGCGATATCG");

	}
}

TEST_CASE("Testing Fastq::nextRecord") {

	Gz::Reader reader("test_data/test.sub.1.fq.gz");
	std::optional<Fastq::Rec> rec = Fastq::nextRecord(reader);
	CHECK(rec->seq_id == "V350082487L3C001R0020000004/1");
	CHECK(rec->seq == "GAACTCTTAGACGGTGCAAGCGCAGAATTTGACATGGATCTTGTATCAAAGGGAGAACTTTCACCTGTATTTTTCGGTTCTGCACTGACAAATTTTGGTGTGGAAACATTTTTAAAGCATTTTTTACAGATGACATATTCTCCATTGCCA");
	CHECK(rec->qual == "6187EBB=F5<?61?A@E</9>BE>E:6ED:8BFF1==DD@D@CEE?FF@=074EBDEFFCCFE7FE8EFEFDB:<DEFE3DDEEFD5EFFFEAFCF5FF7F;/FA&FEA?F=CDCAEFFD@D3FFEFFFF??EFFGFFBFFBFF4DC:D");


	std::vector<std::optional<Fastq::Rec>> recs;

	while (rec) {
		recs.push_back(rec);
		rec = Fastq::nextRecord(reader);
	}

	CHECK(recs.size() == 2500);
	auto last_rec = recs[recs.size()-1];
	CHECK(last_rec->seq_id == "V350082487L3C001R0020157281/1");
	CHECK(last_rec->seq == "GGACATCAGATGATATCCCTCGCCTACGGCGCCAAGACATATAAGCTCAAGTTTGGTCACCGCGGCGGAAACCACCCGGTCATGAACCTCGACACCAACAAGATAGAGATAACCTCGCAGAACCACAGCTATGCCGTTGACCCGAAAACG");
	CHECK(last_rec->qual == "@C,CBA;>E8ECB<<E9?A77>=?F4==@95FAE>,F9DFFF@8EEF?FDA0?E:5D@F)EC<E>A:@,=DEDFCEDADD?D8=DEC@FFD'EFDAFF@FEE;FFC7E:DE9FEBEDFEB>DFDD<99FF>DCDFD<AE=;9ECDEAD84");

	// loops
	Gz::Reader nonexistent_reader("test_data/nonexistent");
	auto empty_rec = Fastq::nextRecord(nonexistent_reader);
	CHECK(!empty_rec);
}

TEST_CASE("Testing Fastq::nextRecordPair") {

	Gz::Reader reader1("test_data/test.sub.1.fq.gz");
	Gz::Reader reader2("test_data/test.sub.2.fq.gz");
	
	std::optional<Fastq::Pair> rec_pair = Fastq::nextRecordPair(reader1, reader2);
	CHECK(rec_pair->first.seq_id == "V350082487L3C001R0020000004/1");
	CHECK(rec_pair->first.seq == "GAACTCTTAGACGGTGCAAGCGCAGAATTTGACATGGATCTTGTATCAAAGGGAGAACTTTCACCTGTATTTTTCGGTTCTGCACTGACAAATTTTGGTGTGGAAACATTTTTAAAGCATTTTTTACAGATGACATATTCTCCATTGCCA");
	CHECK(rec_pair->first.qual == "6187EBB=F5<?61?A@E</9>BE>E:6ED:8BFF1==DD@D@CEE?FF@=074EBDEFFCCFE7FE8EFEFDB:<DEFE3DDEEFD5EFFFEAFCF5FF7F;/FA&FEA?F=CDCAEFFD@D3FFEFFFF??EFFGFFBFFBFF4DC:D");

	CHECK(rec_pair->second.seq_id == "V350082487L3C001R0020000004/2");
	CHECK(rec_pair->second.seq == "GGGTAACTTCCATGCCTGCAGTAAATTTACCGGAACAGATACGCATAAATGCGATCCTGTCCCTGTGGTTTTTATTCATATTTGCCTGAATTTTAAAAACAAATGCGGAAAAATCTTCATCAAACGGATTTTTTTCTCCCTCATTTGACC");
	CHECK(rec_pair->second.qual == ";D<=ECCDFCEEEEBEEF;E@DFFFEECFD8C=EE=EF>EF=C=FEFEDFEDE>E@DFEDBFEFF6FEEEFCCFEFB=DF9EFE<EFEBEEECAACDBFAF=>BD/=7CEFFFB9E5DDF6EAF;/C4@EF@FDFEF4BFE-A>EFD=@C");

	std::vector<std::optional<Fastq::Pair>> rec_pairs;

	while (rec_pair) {
		rec_pairs.push_back(rec_pair);
		rec_pair = Fastq::nextRecordPair(reader1, reader2);
	}

	CHECK(rec_pairs.size() == 2500);
	auto last_rec = rec_pairs[rec_pairs.size()-1];
	CHECK(last_rec->first.seq_id == "V350082487L3C001R0020157281/1");
	CHECK(last_rec->first.seq == "GGACATCAGATGATATCCCTCGCCTACGGCGCCAAGACATATAAGCTCAAGTTTGGTCACCGCGGCGGAAACCACCCGGTCATGAACCTCGACACCAACAAGATAGAGATAACCTCGCAGAACCACAGCTATGCCGTTGACCCGAAAACG");
	CHECK(last_rec->first.qual == "@C,CBA;>E8ECB<<E9?A77>=?F4==@95FAE>,F9DFFF@8EEF?FDA0?E:5D@F)EC<E>A:@,=DEDFCEDADD?D8=DEC@FFD'EFDAFF@FEE;FFC7E:DE9FEBEDFEB>DFDD<99FF>DCDFD<AE=;9ECDEAD84");

	CHECK(last_rec->second.seq_id == "V350082487L3C001R0020157281/2");
	CHECK(last_rec->second.seq == "ACGAGCACCTTTTTTATATCAGTTCTCTTAGGCACGCTCCTTCGCCTCCTTCATATAGTCGATAAACTGACCGAACAGATAGCTGCTGTCCTGCGGACCGGGCGCGCTCTCGGGGTGATACTGCACGCTGAACACAAGGACTTTCTTGCA");
	CHECK(last_rec->second.qual == "FFFAFAFBEFE6DDEFEFF7FE3F@D@AFD@AED>FBFDDFFBF?A0:FFFEFAFEAF90F?EBFFFFCE?EFDE8FFCFAD81ACFGBDBCF:;FE<BBA7BF4F1D3@=CBFF>FCECF/E=D-F2?F-F7F>=F99-FEF9@F=5E=");
}

TEST_CASE("testing hasing") {
	std::string t1 = "nnaCAGCAGTAAAAGCTAAAAGAACGAATACCACaga";
	size_t hash_n = 1;
	size_t kmer_size = 31;
	ntHashIterator itr(t1, hash_n, kmer_size);
	++itr;
	++itr;
	++itr;
	uint64_t hash_value = (*itr)[0];
	//CHECK(hash_value == 1489897631772959496);
	std::string t2 = "AAAcagcagtaaaagctaaaagaacgaataccacAGA";
	ntHashIterator itr2(t2, hash_n, kmer_size);
	++itr2;
	++itr2;
	++itr2;
	uint64_t t2_hash_value = (*itr2)[0];
	CHECK(t2_hash_value == 1489897631772959496);

}
