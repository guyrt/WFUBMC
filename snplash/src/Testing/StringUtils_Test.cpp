#include <gtest/gtest.h>
#include "../engine/utils/stringutils.h"

// Test spaced_string

TEST(StringPadding, HandlesFrontAndBackPadding) {
	std::string s = strnutils::spaced_string("test",0);
	ASSERT_STREQ("test",s.c_str()) << "Test of print to empty string";

	s = strnutils::spaced_string("test",5);
	ASSERT_STREQ("test ",s.c_str()) << "Test of print to extra string";

	s = strnutils::spaced_string("test",4,2);
	ASSERT_STREQ("  test",s.c_str()) << "Test of format string with left padding"; 
}

TEST(IntegerLog, IntLogCorrect) {
	ASSERT_EQ(0,strnutils::int_log(0)) << "Test of log(0)";
	ASSERT_EQ(0,strnutils::int_log(1)) << "Test of log(1)";	
	ASSERT_EQ(0,strnutils::int_log(9)) << "Test of log(9)";
	ASSERT_EQ(1,strnutils::int_log(10)) << "Test of log(10)";
	ASSERT_EQ(1,strnutils::int_log(11)) << "Test of log(11)";

	ASSERT_EQ(0,strnutils::int_log(-1)) << "Test of log(-1)";
	ASSERT_EQ(0,strnutils::int_log(-9)) << "Test of log(-9)";
	ASSERT_EQ(1,strnutils::int_log(-10)) << "Test of log(-10)";
	ASSERT_EQ(1,strnutils::int_log(-11)) << "Test of log(-11)";

	ASSERT_EQ(1,strnutils::int_log(99)) << "Test of log(99)";
	ASSERT_EQ(2,strnutils::int_log(100)) << "Test of log(100)";
	ASSERT_EQ(2,strnutils::int_log(200)) << "Test of log(200)";
	ASSERT_EQ(3,strnutils::int_log(1001)) << "Test of log(1001)";
	ASSERT_EQ(4,strnutils::int_log(10001)) << "Test of log(10001)";
	ASSERT_EQ(5,strnutils::int_log(100001)) << "Test of log(100001)";
}


TEST(IntegerLogDouble, IntLogCorrectForDouble) {
	ASSERT_EQ(0,strnutils::int_log(0.0)) << "Test of log(0.0)";
	ASSERT_EQ(0,strnutils::int_log(1.0)) << "Test of log(1.0)";	
	ASSERT_EQ(0,strnutils::int_log(9.0)) << "Test of log(9.0)";
	ASSERT_EQ(1,strnutils::int_log(10.0)) << "Test of log(10.0)";
	ASSERT_EQ(1,strnutils::int_log(11.0)) << "Test of log(11.0)";

	ASSERT_EQ(0,strnutils::int_log(-1.0)) << "Test of log(-1.0)";
	ASSERT_EQ(0,strnutils::int_log(-9.0)) << "Test of log(-9.0)";
	ASSERT_EQ(1,strnutils::int_log(-10.0)) << "Test of log(-10.0)";
	ASSERT_EQ(1,strnutils::int_log(-11.0)) << "Test of log(-11.0)";

	ASSERT_EQ(1,strnutils::int_log(99.0)) << "Test of log(99.0)";
	ASSERT_EQ(2,strnutils::int_log(100.0)) << "Test of log(100.0)";
	ASSERT_EQ(2,strnutils::int_log(200.0)) << "Test of log(200.0)";
	ASSERT_EQ(3,strnutils::int_log(1001.0)) << "Test of log(1001.0)";
	ASSERT_EQ(4,strnutils::int_log(10001.0)) << "Test of log(10001.0)";
	ASSERT_EQ(5,strnutils::int_log(100001.0)) << "Test of log(100001.0)";
}
