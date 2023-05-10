#include "doctest.h"
#include "utils.h"
#include <cstdio>

TEST_CASE("Test Gz::trimNewlineInplace") {
	std::string t1 = "";
	trimNewlineInplace(t1);
	CHECK(t1 == "");

	std::string t2 = "\n";
	trimNewlineInplace(t2);
	CHECK(t2 == "");

	std::string t3 = "\r\n";
	trimNewlineInplace(t3);
	CHECK(t3 == "");

	std::string t4 = "atg";
	trimNewlineInplace(t4);
	CHECK(t4 == "atg");

	std::string t5 = "atg\n";
	trimNewlineInplace(t5);
	CHECK(t5 == "atg");

	std::string t6 = "atg\r\n";
	trimNewlineInplace(t6);
	CHECK(t6 == "atg");
}

TEST_CASE("Test Gz::Writer::writeLine") {
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
