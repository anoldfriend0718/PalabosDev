#include "../src/parameters/singlePhaseFlowParameters.h"
#include "core/globalDefs.h"
#include "palabos2D.h"
#include "palabos2D.hh"
#include "gtest/gtest.h"

using namespace testing;
using namespace plb;

typedef double T;
class SinglePhaseFlowParameterExample {
public:
  SinglePhaseFlowParameterExample() {
    T physicalU = 0.5;
    pluint latticeCharacteristicLength = 100;
    T resolution = 1e-6;
    T re = 5;
    T latticeU = 0.01;
    pluint latticeLx = 100;
    pluint latticeLy = 100;
    flowParam = SinglePhaseFlowParam<T>(physicalU, latticeU, re,
                                        latticeCharacteristicLength, resolution,
                                        latticeLx, latticeLy);
  }
  SinglePhaseFlowParam<T> flowParam;
};

class SinglePhaseFlowParameterTest : public Test {
public:
  SinglePhaseFlowParam<T> flowParam =
      SinglePhaseFlowParameterExample().flowParam;
};

TEST_F(SinglePhaseFlowParameterTest, TestIfGetCorrectParameters) {
  ASSERT_DOUBLE_EQ(flowParam.getPhysicalCharacteristicLength(), 1e-4);
  ASSERT_DOUBLE_EQ(flowParam.getPhysicalLx(), 1e-4);
  ASSERT_DOUBLE_EQ(flowParam.getPhysicalLy(), 1e-4);
  ASSERT_DOUBLE_EQ(flowParam.getDeltaX(), 1e-6);
  ASSERT_DOUBLE_EQ(flowParam.getDeltaT(), 2.0 * 1e-8);
  ASSERT_DOUBLE_EQ(flowParam.getLatticeNu(), 0.2);
  ASSERT_DOUBLE_EQ(flowParam.getTau(), 1.1);
  ASSERT_DOUBLE_EQ(flowParam.getOmega(), 1.0 / 1.1);
}
