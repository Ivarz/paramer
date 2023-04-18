#include "doctest.h"
#include "utils.h"

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
