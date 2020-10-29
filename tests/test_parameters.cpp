#include "core/globalDefs.h"
#include "palabos2D.h"
#include "palabos2D.hh"
#include "gtest/gtest.h"

using namespace testing;
using namespace plb;

/// Numeric parameters for isothermal, incompressible flow.
template <typename T> class SinglePhaseFlowParam {
public:
  SinglePhaseFlowParam(T physicalU_, T latticeU_, T Re_,
                       pluint latticeCharacteristicLength_,
                       T physicalLengthPerPixel, pluint latticeLx_,
                       pluint LatticeLy_, pluint latticeLz_ = T())
      : physicalU(physicalU_), latticeU(latticeU_), Re(Re_),
        latticeCharacteristicLength(latticeCharacteristicLength_),
        resolution(physicalLengthPerPixel), latticeLx(latticeLx_),
        latticeLy(LatticeLy_), latticeLz(latticeLz_) {}
  SinglePhaseFlowParam() = default;

  /// velocity in lattice units (proportional to Mach number)
  T getLatticeU() const { return latticeU; }
  /// velocity in physical units
  T getPhysicalU() const { return physicalU; }
  /// Reynolds number
  T getRe() const { return Re; }
  /// physical characteristic length
  T getPhysicalCharacteristicLength() const {
    return latticeCharacteristicLength * resolution;
  }
  /// resolution
  T getResolution() const { return resolution; }
  pluint getLatticeCharacteristicLength() const {
    return latticeCharacteristicLength;
  }
  /// x-length in physical units
  T getPhysicalLx() const { return resolution * latticeLx; }
  /// y-length in physical units
  T getPhysicalLy() const { return resolution * latticeLy; }
  /// z-length in physical units
  T getPhysicalLz() const { return resolution * latticeLz; }
  /// lattice spacing in dimensionless units
  T getDeltaX() const { return resolution; }
  /// time step in dimensionless units
  T getDeltaT() const { return getDeltaX() * getLatticeU() / getPhysicalU(); }
  /// conversion from dimensionless to lattice units for space coordinate
  pluint nCell(T l) const { return (int)(l / getDeltaX() + (T)0.5); }
  /// conversion from dimensionless to lattice units for time coordinate
  pluint nStep(T t) const { return (int)(t / getDeltaT() + (T)0.5); }
  /// number of lattice cells in x-direction
  pluint getNx(bool offLattice = false) const {
    return nCell(getPhysicalLx()) + 1 + (int)offLattice;
  }
  /// number of lattice cells in y-direction
  pluint getNy(bool offLattice = false) const {
    return nCell(getPhysicalLy()) + 1 + (int)offLattice;
  }
  /// number of lattice cells in z-direction
  pluint getNz(bool offLattice = false) const {
    return nCell(getPhysicalLz()) + 1 + (int)offLattice;
  }
  /// viscosity in lattice units
  T getLatticeNu() const {
    return getLatticeU() * (T)getLatticeCharacteristicLength() / Re;
  }
  /// relaxation time
  T getTau() const { return (T)3 * getLatticeNu() + (T)0.5; }
  /// relaxation frequency
  T getOmega() const { return (T)1 / getTau(); }

private:
  T physicalU, latticeU, Re;
  pluint latticeCharacteristicLength;
  T resolution;
  pluint latticeLx, latticeLy, latticeLz;
};

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
