#include "palabos2D.h"
#include "palabos2D.hh"
#include "gtest/gtest.h"

using namespace testing;
using namespace plb;

int main(int argc, char **argv) {
  plbInit(&argc, &argv); // should init plb first, otherwise undefined behavior
                         // happen when using the Palabos data structure
  InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
