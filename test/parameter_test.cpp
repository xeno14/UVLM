#include <gtest/gtest.h>

#include "parameter.h"

DEFINE_PARAM_int(int_default, 100);
TEST(DefaultValueTest, int) {
  EXPECT_EQ(100, PARAM_int_defalut); 
}
