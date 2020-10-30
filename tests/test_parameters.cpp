#include "../src/parameters/singlePhaseFlowParameters.h"
#include "basicDynamics/isoThermalDynamics.h"
#include "core/dynamicsIdentifiers.h"
#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "palabos2D.h"
#include "palabos2D.hh"
#include "gtest/gtest.h"
#include <plog/Severity.h>

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

enum DynamicsName { BGK_Ma2 = 0, RM, HM, CM, CHM, K, GH, RR };

enum HOOmega { SRT = 0, REG };

template <typename T, template <typename U> class Descriptor> class LBMModel2D {
public:
  LBMModel2D(DynamicsName dynamics, HOOmega hoOmega, T omega_) : omega(omega_) {
    setDynamics(dynamics);
    setOmega(hoOmega);
  }

  std::unique_ptr<Dynamics<T, Descriptor>> getDynamics() const {
    std::unique_ptr<Dynamics<T, Descriptor>> dynamics(
        isoThermalDynamics->clone());
    return dynamics;
  }

  Array<T, Descriptor<T>::numRelaxationTimes> getAllOmega() const {
    return allOmega;
  }

private:
  std::unique_ptr<IsoThermalBulkDynamics<T, Descriptor>> isoThermalDynamics;
  T omega;
  Array<T, Descriptor<T>::numRelaxationTimes> allOmega;

  void setDynamics(DynamicsName &dynamics) {
    switch (dynamics) {
    case (BGK_Ma2):
      isoThermalDynamics.reset(new BGKdynamics<T, Descriptor>(omega));
      break;
    case (RM):
      isoThermalDynamics.reset(new RMdynamics<T, Descriptor>(omega));
      break;
    case (HM):
      isoThermalDynamics.reset(new HMdynamics<T, Descriptor>(omega));
      break;
    case (CM):
      isoThermalDynamics.reset(new CMdynamics<T, Descriptor>(omega));
      break;
    case (CHM):
      isoThermalDynamics.reset(new CHMdynamics<T, Descriptor>(omega));
      break;
    case (GH):
      isoThermalDynamics.reset(new GHdynamics<T, Descriptor>(omega));
      break;
    case (RR):
      isoThermalDynamics.reset(new RRdynamics<T, Descriptor>(omega));
      break;
    case (K):
      isoThermalDynamics.reset(new Kdynamics<T, Descriptor>(omega));
      break;
    default:
      PLOG(plog::error) << "Error: dynamics name does not exist." << std::endl;
      exit(-1);
      break;
    }
  }
  void setOmega(HOOmega hoOmega) {
    PLB_PRECONDITION(Descriptor<T>::numRelaxationTimes == 4);
    switch (hoOmega) {
    case (SRT):
      allOmega[0] = omega;
      allOmega[1] = omega;
      allOmega[2] = omega;
      allOmega[3] = omega;
      break;
    case (REG):
      allOmega[0] = omega;
      allOmega[1] = omega;
      allOmega[2] = 1;
      allOmega[3] = 1;
      break;
    default:
      PLOG(plog::error) << "Error: high order omega catelog does not exist."
                        << std::endl;
      exit(-1);
      break;
    }
  }
};

#define DESCRIPTOR descriptors::D2Q9Descriptor

template <typename T_, template <typename U> class Descriptor>
std::string getDynamicsName(LBMModel2D<T_, Descriptor> &model) {
  plb::meta::DynamicsRegistration<double, DESCRIPTOR> &dynamicsRegistration =
      plb::meta::dynamicsRegistration<T_, DESCRIPTOR>();

  auto id = model.getDynamics()->getId();
  std::string dynamicsName = dynamicsRegistration.getName(id);
  return dynamicsName;
};

TEST(LBMModelParameter, TestIfCorrectDynamicsCatelog) {

  LBMModel2D<T, DESCRIPTOR> lbmModel1(BGK_Ma2, SRT, 1.0);
  ASSERT_EQ(getDynamicsName(lbmModel1), "BGK");

  LBMModel2D<T, DESCRIPTOR> lbmModel2(RM, SRT, 1.0);
  ASSERT_EQ(getDynamicsName(lbmModel2), "RM");

  LBMModel2D<T, DESCRIPTOR> lbmModel3(HM, SRT, 1.0);
  ASSERT_EQ(getDynamicsName(lbmModel3), "HM");

  LBMModel2D<T, DESCRIPTOR> lbmModel4(CM, SRT, 1.0);
  ASSERT_EQ(getDynamicsName(lbmModel4), "CM");

  LBMModel2D<T, DESCRIPTOR> lbmModel5(CHM, SRT, 1.0);
  ASSERT_EQ(getDynamicsName(lbmModel5), "CHM");

  LBMModel2D<T, DESCRIPTOR> lbmModel6(K, SRT, 1.0);
  ASSERT_EQ(getDynamicsName(lbmModel6), "K");

  LBMModel2D<T, DESCRIPTOR> lbmModel7(GH, SRT, 1.0);
  ASSERT_EQ(getDynamicsName(lbmModel7), "GH");

  LBMModel2D<T, DESCRIPTOR> lbmModel8(RR, SRT, 1.0);
  ASSERT_EQ(getDynamicsName(lbmModel8), "RR");
}

TEST(LBMModelParameter,TestIfCorrectAllOmega)
{
  LBMModel2D<T, DESCRIPTOR> lbmModel1(CHM, SRT, 1.2);
  ASSERT_EQ(lbmModel1.getAllOmega()[0], 1.2);
  ASSERT_EQ(lbmModel1.getAllOmega()[1], 1.2);
  ASSERT_EQ(lbmModel1.getAllOmega()[2], 1.2);
  ASSERT_EQ(lbmModel1.getAllOmega()[3], 1.2);

  LBMModel2D<T, DESCRIPTOR> lbmModel2(CHM, REG, 1.2);
  ASSERT_EQ(lbmModel2.getAllOmega()[0], 1.2);
  ASSERT_EQ(lbmModel2.getAllOmega()[1], 1.2);
  ASSERT_EQ(lbmModel2.getAllOmega()[2], 1.0);
  ASSERT_EQ(lbmModel2.getAllOmega()[3], 1.0);
}
// CHM, K, GH, RR