#include "doctest.h"
#include "kraken2.h"

TEST_CASE("Test Kraken2::nextRecord") {
	{
		Gz::Reader gzreader("test_data/empty_k2.out.txt");
		std::optional<Kraken2::Rec> rec1 = Kraken2::nextRecord(gzreader);
		CHECK(!rec1);
	}
	{
		Gz::Reader gzreader("test_data/t1.k2out.txt.gz");
		std::optional<Kraken2::Rec> rec1 = Kraken2::nextRecord(gzreader);
		CHECK( rec1->seq_id == "EVEC_scaffold0000001");
		CHECK( rec1->taxid == "0");
		CHECK( !rec1->paired_end);
	}
	{
		Gz::Reader gzreader("test_data/t2.k2out.txt");
		std::optional<Kraken2::Rec> rec1 = Kraken2::nextRecord(gzreader);
		CHECK( rec1->seq_id == "EVEC_scaffold0000001");
		CHECK( rec1->taxid == "0");
		CHECK( !rec1->paired_end);
	}
	{
		Gz::Reader gzreader("test_data/t3.k2out.txt");
		std::optional<Kraken2::Rec> rec1 = Kraken2::nextRecord(gzreader);
		CHECK( rec1->seq_id == "EVEC_scaffold0000001");
		CHECK( rec1->taxid == "0");
		CHECK( !rec1->paired_end);
		std::vector<std::string> true_taxids = {"3", "0", "A", "0", "4", "0", "A", "0"};
		std::vector<size_t> true_kmers = {16, 327, 160, 1236, 16, 100, 126, 5089};
		auto r1_kmers = rec1->getR1Kmers();
		for (size_t i=0; i < true_taxids.size(); i++) {
			CHECK(true_taxids[i] == r1_kmers[i].first);
			CHECK(true_kmers[i] == r1_kmers[i].second);
		}

	}
}
