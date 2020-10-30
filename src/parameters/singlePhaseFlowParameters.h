#include "core/globalDefs.h"
#include "io/parallelIO.h"
#include "plog/Log.h"
#include "plog/Severity.h"
#include <iostream>
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
  /// relaxation time for D2Q9 with cs*cs=1/3
  T getTau() const { return (T)3 * getLatticeNu() + (T)0.5; }
  /// relaxation frequency
  T getOmega() const { return (T)1 / getTau(); }

  bool writeLogFile(std::string const &title = "") {
    std::string fullName =
        global::directories().getLogOutDir() + "singlePhaseFlowParameters.dat";
    plb_ofstream ofile(fullName.c_str());
    if (!ofile.is_open()) {
      pcerr << "could not open the log file: " << fullName << std::endl;
      return false;
    }
    ofile << title << "\n\n";
    ofile << "Reynolds number:                    Re=" << getRe() << "\n";
    ofile << "number of lattice cells in x-axis   nx=" << getNx() << "\n";
    ofile << "number of lattice cells in y-axis   ny=" << getNy() << "\n";
    ofile << "number of lattice cells in z-axis   nz=" << getNz() << "\n";
    ofile << "Velocity in lattice units:          u=" << getLatticeU() << "\n";
    ofile << "Relaxation frequency:               omega=" << getOmega() << "\n";
    ofile << "Relaxation time:                    tau=" << getTau() << "\n";
    ofile << "Extent of the physical system:      lx=" << getPhysicalLx()
          << "\n";
    ofile << "Extent of the physical system:      ly=" << getPhysicalLy()
          << "\n";
    ofile << "Extent of the physical system:      lz=" << getPhysicalLz()
          << "\n";
    ofile << "Grid spacing deltaX:                dx=" << getDeltaX() << "\n";
    ofile << "Time step deltaT:                   dt=" << getDeltaT() << "\n";
    PLOG(plog::info) << "single phase flow parameter is output to " << fullName;
    return true;
  }

private:
  T physicalU, latticeU, Re;
  pluint latticeCharacteristicLength;
  T resolution;
  pluint latticeLx, latticeLy, latticeLz;
};
