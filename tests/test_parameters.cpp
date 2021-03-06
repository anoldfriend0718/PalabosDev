
#include "palabos2D.h"
#include "palabos2D.hh"
#include "gtest/gtest.h"
#include <plog/Severity.h>
#include "parameters/singlePhaseFlowParameters.h"
#include "parameters/LBMModelParser2D.h"
#include "parameters/LBMModelParser2D.hh"
#include "basicDynamics/isoThermalDynamics.h"
#include "basicDynamics/comprehensiveIsoThermalDynamics.h"
#include "customizedutil/utilities.h"

using namespace testing;
using namespace plb;

#define DESCRIPTOR descriptors::D2Q9Descriptor

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



TEST(LBMModelParameter, TestIfCorrectDynamicsCatelog) {

  LBMModelParser2D<T, DESCRIPTOR> lbmModel1(BGK_Ma2, SRT, 1.0);
  ASSERT_EQ(util::getDynamicsName(lbmModel1), "BGK");

  LBMModelParser2D<T, DESCRIPTOR> lbmModel2(RM, SRT, 1.0);
  ASSERT_EQ(util::getDynamicsName(lbmModel2), "RM");

  LBMModelParser2D<T, DESCRIPTOR> lbmModel3(HM, SRT, 1.0);
  ASSERT_EQ(util::getDynamicsName(lbmModel3), "HM");

  LBMModelParser2D<T, DESCRIPTOR> lbmModel4(CM, SRT, 1.0);
  ASSERT_EQ(util::getDynamicsName(lbmModel4), "CM");

  LBMModelParser2D<T, DESCRIPTOR> lbmModel5(CHM, SRT, 1.0);
  ASSERT_EQ(util::getDynamicsName(lbmModel5), "CHM");

  LBMModelParser2D<T, DESCRIPTOR> lbmModel6(K, SRT, 1.0);
  ASSERT_EQ(util::getDynamicsName(lbmModel6), "K");

  LBMModelParser2D<T, DESCRIPTOR> lbmModel7(GH, SRT, 1.0);
  ASSERT_EQ(util::getDynamicsName(lbmModel7), "GH");

  LBMModelParser2D<T, DESCRIPTOR> lbmModel8(RR, SRT, 1.0);
  ASSERT_EQ(util::getDynamicsName(lbmModel8), "RR");
}

TEST(LBMModelParameter,TestIfGetCorrectAllOmega)
{
  LBMModelParser2D<T, DESCRIPTOR> lbmModel1(CHM, SRT, 1.2);
  ASSERT_EQ(lbmModel1.getAllOmega()[0], 1.2);
  ASSERT_EQ(lbmModel1.getAllOmega()[1], 1.2);
  ASSERT_EQ(lbmModel1.getAllOmega()[2], 1.2);
  ASSERT_EQ(lbmModel1.getAllOmega()[3], 1.2);

  LBMModelParser2D<T, DESCRIPTOR> lbmModel2(CHM, REG, 1.2);
  ASSERT_EQ(lbmModel2.getAllOmega()[0], 1.2);
  ASSERT_EQ(lbmModel2.getAllOmega()[1], 1.2);
  ASSERT_EQ(lbmModel2.getAllOmega()[2], 1.0);
  ASSERT_EQ(lbmModel2.getAllOmega()[3], 1.0);
}

TEST(LBMModelParameter,TestIfSetCorrectAllOmega)
{
  LBMModelParser2D<T, DESCRIPTOR> lbmModel2(CHM, REG, 1.2);
  plb::Array<T,DESCRIPTOR<T>::numRelaxationTimes> allOmega = CHMdynamics<T, DESCRIPTOR>::allOmega;
  ASSERT_EQ(allOmega[0], 1.2);
  ASSERT_EQ(allOmega[1], 1.2);
  ASSERT_EQ(allOmega[2], 1.0);
  ASSERT_EQ(allOmega[3], 1.0);
}

TEST(LBMModelParameter, TestIfCorrectConstructionByString)
{
  LBMModelParser2D<T, DESCRIPTOR> lbmModel1("BGK_Ma2","SRT",1.2);
  ASSERT_EQ(util::getDynamicsName(lbmModel1), "BGK");
  ASSERT_EQ(lbmModel1.getAllOmega()[0], 1.2);

  LBMModelParser2D<T, DESCRIPTOR> lbmModel2("RR","REG",1.2);
  ASSERT_EQ(util::getDynamicsName(lbmModel2), "RR");
  plb::Array<T, DESCRIPTOR<T>::numRelaxationTimes> allOmega =
      RRdynamics<T, DESCRIPTOR>::allOmega;
  ASSERT_EQ(allOmega[0], 1.2);
  ASSERT_EQ(allOmega[1], 1.2);
  ASSERT_EQ(allOmega[2], 1.0);
  ASSERT_EQ(allOmega[3], 1.0);


}
