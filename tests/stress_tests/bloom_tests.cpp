#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "bloom.h"
#include <cstdio>

const size_t TARGET_SIZE = 10000000000;
TEST_CASE("Test write and load") {
	{{
		CHECK(true);
		Bloom::Filter t0_bloom = Bloom::Filter(TARGET_SIZE, 31, 1);
		t0_bloom.write("test_data/t0.stress.blm");
	}
	{
		std::optional<Bloom::Filter> t0_bloom_load = Bloom::Filter::load("test_data/t0.stress.blm");
		std::remove("test_data/t0.stress.blm");
		CHECK(t0_bloom_load);
		CHECK(t0_bloom_load->size() == TARGET_SIZE);
		CHECK(t0_bloom_load->setBitsCount() == 0);
		CHECK(t0_bloom_load->kmerSize() == 31);
		CHECK(t0_bloom_load->hashN() == 1);
	}}

	{
		 std::vector<uint64_t> rand_idx = {
			 506185150 , 2070281628 , 1647099284 , 1017519544 , 2893929445 , 4361684722 , 5140858614 , 8536299639 , 7186193590 , 4237109865 ,
			 3553449885 , 8058684656 , 3128432620 , 5589428357 , 9504155698 , 2017373019 , 7084215493 , 3316070372 , 1604151810 , 4452380177 ,
			 7346241811 , 5438887893 , 2883249939 , 1167510616 , 8666193866 , 9959088242 , 5821371180 , 5235751719 , 6826639070 , 7098744644 ,
			 7081094882 , 7058101389 , 4384718948 , 5458652404 , 4624501647 , 2680710576 , 5832874146 , 2159786045 , 3563285325 , 3661031049 ,
			 2698209628 , 3693771946 , 2153893332 , 976426583 , 6458330187 , 4141248158 , 3860548612 , 3251600135 , 8996508328 , 6896058760 ,
			 3462473717 , 5686132308 , 9903449756 , 6924015291 , 3151119323 , 2084270746 , 8868738517 , 1858407557 , 4528115593 , 3108924456 ,
			 3182279604 , 2028256608 , 2547032007 , 1209188961 , 2153083547 , 6715527497 , 8453509042 , 4579946273 , 2440143205 , 6317044506 ,
			 940305672 , 2573875159 , 3365128928 , 9281512511 , 4537281737 , 374008263 , 3372773976 , 4586379048 , 3185535620 , 9323476936 ,
			 8191735 , 8997014354 , 3315094186 , 9452081852 , 8603695989 , 7558197788 , 6298709369 , 5990563177 , 4320300068 , 6482917669 ,
			 903028172 , 6631252666 , 8489607383 , 8122639148 , 19829609 , 1258499368 , 1674532180 , 6482000094 , 5124930475 , 6481629739
		 };
		std::vector<uint8_t> rand_bytes = {
			109 , 244 , 251 , 11 , 131 , 198 , 158 , 242 , 185 , 180 ,
			16 , 208 , 75 , 152 , 43 , 121 , 139 , 208 , 80 , 49 ,
			212 , 246 , 194 , 163 , 165 , 102 , 187 , 174 , 12 , 200 ,
			150 , 53 , 124 , 16 , 144 , 99 , 213 , 177 , 72 , 182 ,
			239 , 91 , 157 , 170 , 119 , 234 , 206 , 5 , 84 , 7 ,
			242 , 227 , 141 , 196 , 205 , 29 , 255 , 44 , 117 , 227 ,
			200 , 152 , 45 , 108 , 48 , 53 , 237 , 187 , 138 , 167 ,
			52 , 10 , 209 , 99 , 11 , 177 , 11 , 228 , 134 , 76 ,
			236 , 142 , 223 , 187 , 131 , 195 , 41 , 184 , 114 , 108 ,
			50 , 244 , 3 , 143 , 41 , 180 , 238 , 247 , 130 , 210 
		};
		{

			std::optional<Bloom::Filter> t1_bloom = Bloom::Filter(TARGET_SIZE, 31, 1);
			for (size_t i=0; i < rand_idx.size(); i++) {
				size_t idx = rand_idx[i];
				size_t value = rand_bytes[i];
				t1_bloom->at(idx) = value;
			}
			t1_bloom->at(0) = 1;
			t1_bloom->at(3) = 3;
			t1_bloom->at(8) = 2;
			t1_bloom->at(10) = 4;
			t1_bloom->at(TARGET_SIZE-1) = 255;
			t1_bloom->write("test_data/t1.stress.blm");
		}
		{
			std::optional<Bloom::Filter> t1_bloom_load = Bloom::Filter::load("test_data/t1.stress.blm");
			std::remove("test_data/t1.stress.blm");
			CHECK(t1_bloom_load);
			CHECK(t1_bloom_load->at(0) == 1);
			CHECK(t1_bloom_load->at(3) == 3);
			CHECK(t1_bloom_load->at(8) == 2);
			CHECK(t1_bloom_load->at(10) == 4);
			CHECK(t1_bloom_load->at(TARGET_SIZE-1) == 255);
			for (size_t i=0; i < rand_idx.size(); i++) {
				size_t idx = rand_idx[i];
				size_t value = rand_bytes[i];
				CHECK(t1_bloom_load->at(idx) == value);
			}
			std::cerr << "False positive rate: " << t1_bloom_load->falsePostiveRate() << '\n';
		}
	}
}

TEST_CASE("Test write and load raw") {
	{
	Bloom::Filter t0_bloom = Bloom::Filter(TARGET_SIZE, 31, 1);
	t0_bloom.writeRaw("test_data/t0.stress.blm");
	}
	{
	std::optional<Bloom::Filter> t0_bloom_load = Bloom::Filter::loadRaw("test_data/t0.stress.blm");
	std::remove("test_data/t0.stress.blm");
	CHECK(t0_bloom_load->size() == TARGET_SIZE);
	CHECK(t0_bloom_load->setBitsCount() == 0);
	CHECK(t0_bloom_load->kmerSize() == 31);
	CHECK(t0_bloom_load->hashN() == 1);
	}

	{
	std::optional<Bloom::Filter> t1_bloom = Bloom::Filter(TARGET_SIZE, 31, 1);
	t1_bloom->at(0) = 1;
	t1_bloom->at(3) = 3;
	t1_bloom->at(8) = 2;
	t1_bloom->at(10) = 4;
	t1_bloom->at(TARGET_SIZE-1) = 255;
	t1_bloom->writeRaw("test_data/t1.stress.blm");
	}
	{
	std::optional<Bloom::Filter> t1_bloom_load = Bloom::Filter::loadRaw("test_data/t1.stress.blm");
	std::remove("test_data/t1.stress.blm");
	CHECK(t1_bloom_load->at(0) == 1);
	CHECK(t1_bloom_load->at(3) == 3);
	CHECK(t1_bloom_load->at(8) == 2);
	CHECK(t1_bloom_load->at(10) == 4);
	CHECK(t1_bloom_load->at(TARGET_SIZE-1) == 255);
	}
}
