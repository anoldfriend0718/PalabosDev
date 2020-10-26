#include "gtest/gtest.h"
using namespace testing;

int main(int argc, char *argv[]) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

TEST(try, temp) {
  int a = 1;
  ASSERT_EQ(a, 1);
}