#include "doctest.h"
#include "utils.h"
#include <cstdio>

TEST_CASE("Test Utils::trimNewlineInplace") {
	std::string t1 = "";
	Utils::trimNewlineInplace(t1);
	CHECK(t1 == "");

	std::string t2 = "\n";
	Utils::trimNewlineInplace(t2);
	CHECK(t2 == "");

	std::string t3 = "\r\n";
	Utils::trimNewlineInplace(t3);
	CHECK(t3 == "");

	std::string t4 = "atg";
	Utils::trimNewlineInplace(t4);
	CHECK(t4 == "atg");

	std::string t5 = "atg\n";
	Utils::trimNewlineInplace(t5);
	CHECK(t5 == "atg");

	std::string t6 = "atg\r\n";
	Utils::trimNewlineInplace(t6);
	CHECK(t6 == "atg");
}

TEST_CASE("Test Utils::dataSizeToBytes") {
	CHECK(Utils::dataSizeToBytes("") == 0);
	CHECK(Utils::dataSizeToBytes("1000") == 1000);
	CHECK(Utils::dataSizeToBytes("1K") == 1024);
	CHECK(Utils::dataSizeToBytes("1M") == 1024*1024);
	CHECK(Utils::dataSizeToBytes("1G") == 1024*1024*1024);
	CHECK(Utils::dataSizeToBytes("13K") == 13*1024);
	CHECK(Utils::dataSizeToBytes("13M") == 13*1024*1024);
	CHECK(Utils::dataSizeToBytes("13G") == (size_t)13*1024*1024*1024);
}

TEST_CASE("Test Gz::Writer::writeLine") {
	{
		std::string fname = "test_data/gz_writer_test.txt.gz";
		std::string t1 = "test string1";
		std::string t2 = "test string2";
		std::string t3 = "test string3";
		{
			Gz::Writer gzw = Gz::Writer(fname);
			gzw.writeLine(t1);
			gzw.writeLine(t2);
			gzw.writeLine(t3);
		}
		Gz::Reader gzr = Gz::Reader(fname);
		std::string l1 = gzr.nextLine();
		std::string l2 = gzr.nextLine();
		std::string l3 = gzr.nextLine();
		CHECK(t1 == l1);
		CHECK(t2 == l2);
		CHECK(t3 == l3);
		std::remove(fname.c_str());
	}

	{
		std::string fname = "test_data/gz_writer_test.txt.gz";
		std::string t1 = "test string1";
		std::string t2 = "test string2";
		std::string t3 = "test string3";
		{
			std::optional<Gz::Writer> gzw = {};
			gzw.emplace(fname);
			gzw->writeLine(t1);
			gzw->writeLine(t2);
			gzw->writeLine(t3);
		}
		Gz::Reader gzr = Gz::Reader(fname);
		std::string l1 = gzr.nextLine();
		std::string l2 = gzr.nextLine();
		std::string l3 = gzr.nextLine();
		CHECK(t1 == l1);
		CHECK(t2 == l2);
		CHECK(t3 == l3);
		std::remove(fname.c_str());
	}
}
