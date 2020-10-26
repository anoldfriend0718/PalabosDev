
/** \file
 * Flow around a 2D cylinder inside a channel with Guo's parameter, with the
 * creation of a von Karman vortex street. This example makes use of bounce-back
 * nodes to describe the shape of the cylinder. The outlet is modeled through a
 * Neumann (zero velocity-gradient) condition.
 */
#include "customizedutil/residualTracer.h"
#include "customizedutil/residualTracer.hh"
#include "palabos2D.h"
#include "palabos2D.hh"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
#define DESCRIPTOR D2Q9Descriptor

/// Linearly decreasing pressure profile
T poiseuillePressure(plint iX, IncomprFlowParam<T> const &parameters) {
  T Lx = parameters.getNx() - 1;
  T Ly = parameters.getNy() - 1;
  return 8. * parameters.getLatticeNu() * parameters.getLatticeU() / (Ly * Ly) *
         (Lx - (T)iX);
}

/// Convert pressure to density according to ideal gas law
T poiseuilleDensity(plint iX, IncomprFlowParam<T> const &parameters) {
  return poiseuillePressure(iX, parameters) * DESCRIPTOR<T>::invCs2 + (T)1;
}

/// A functional, used to initialize a pressure boundary to constant density
template <typename T> class ConstantDensity {
public:
  ConstantDensity(T density_) : density(density_) {}
  T operator()(plint iX, plint iY) const { return density; }

private:
  T density;
};

template <typename T> class ConstantVelocityAndPoiseuilleDensity {
public:
  ConstantVelocityAndPoiseuilleDensity(IncomprFlowParam<T> parameters_)
      : parameters(parameters_) {}
  void operator()(plint iX, plint iY, T &rho, Array<T, 2> &u) const {
    rho = poiseuilleDensity(iX, parameters);
    u[0] = parameters.getLatticeU();
    u[1] = T();
  }

private:
  IncomprFlowParam<T> parameters;
};
/// A functional, used to instantiate bounce-back nodes at the locations of
/// the cylinder
template <typename T>
class CylinderShapeDomain2D : public plb::DomainFunctional2D {
public:
  CylinderShapeDomain2D(plb::plint cx_, plb::plint cy_, plb::plint radius)
      : cx(cx_), cy(cy_), radiusSqr(plb::util::sqr(radius)) {}
  virtual bool operator()(plb::plint iX, plb::plint iY) const {
    return plb::util::sqr(iX - cx) + plb::util::sqr(iY - cy) <= radiusSqr;
  }
  virtual CylinderShapeDomain2D<T> *clone() const {
    return new CylinderShapeDomain2D<T>(*this);
  }

private:
  plb::plint cx;
  plb::plint cy;
  plb::plint radiusSqr;
};

void cylinderSetup(
    MultiBlockLattice2D<T, DESCRIPTOR> &lattice,
    IncomprFlowParam<T> const &parameters,
    OnLatticeBoundaryCondition2D<T, DESCRIPTOR> &boundaryCondition, plint cx,
    plint cy, plint radius) {
  const plint nx = parameters.getNx();
  const plint ny = parameters.getNy();

  // // create the periodicity in the upper and bottom boundary
  // lattice.periodicity().toggle(1, true);

  // Create Velocity boundary conditions: inlet
  Box2D inlet(0, 0, 0, ny - 1);
  boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, inlet);
  // .. except on right boundary, where we prefer an outflow condition
  //    (zero velocity-gradient).
  Box2D outlet(nx - 1, nx - 1, 0, ny - 1);
  boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, outlet,
                                                          boundary::outflow);
  Box2D top(1, nx - 2, ny - 1, ny - 1);
  boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, top,
                                                          boundary::freeslip);
  Box2D bottom(1, nx - 2, 0, 0);
  boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, bottom,
                                                          boundary::freeslip);

  const Array<T, 2> inletVelocity = {parameters.getLatticeU(), T()};
  setBoundaryVelocity(lattice, lattice.getBoundingBox(), inletVelocity);
  setBoundaryDensity(lattice, outlet, ConstantDensity<T>(1.));
  initializeAtEquilibrium(lattice, lattice.getBoundingBox(),
                          ConstantVelocityAndPoiseuilleDensity<T>(parameters));

  defineDynamics(lattice, lattice.getBoundingBox(),
                 new CylinderShapeDomain2D<T>(cx, cy, radius),
                 new plb::BounceBack<T, DESCRIPTOR>);

  lattice.initialize();
}

void writeGif(MultiBlockLattice2D<T, DESCRIPTOR> &lattice, plint iter) {
  ImageWriter<T> imageWriter("leeloo");
  imageWriter.writeScaledGif(createFileName("u", iter, 6),
                             *computeVelocityNorm(lattice));
  imageWriter.writeScaledGif(createFileName("rho", iter, 6),
                             *computeDensity(lattice));
}

void writeVTK(MultiBlockLattice2D<T, DESCRIPTOR> &lattice,
              IncomprFlowParam<T> const &parameters, plint iter) {
  T dx = parameters.getDeltaX();
  T dt = parameters.getDeltaT();
  VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), dx);
  vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm",
                          dx / dt);
  vtkOut.writeData<2, float>(*computeVelocity(lattice), "velocity", dx / dt);
}

void writeField(const string outputDir,
                MultiBlockLattice2D<T, DESCRIPTOR> &lattice,
                IncomprFlowParam<T> const &parameters, plint iter,
                plint maxdigit = 6) {
  T dx = parameters.getDeltaX();
  T dt = parameters.getDeltaT();
  int precision = 6;
  // density
  const string densityFileName =
      outputDir + createFileName("density", iter, 6) + ".dat";
  plb_ofstream densityFileStream(densityFileName.c_str());
  densityFileStream << setprecision(precision)
                    << *computeDensity(lattice, lattice.getBoundingBox());
  // velocity norm
  const string velocityNormFileName =
      outputDir + createFileName("velocityNorm", iter, maxdigit) + ".dat";
  plb_ofstream velocityNormFileStream(velocityNormFileName.c_str());
  velocityNormFileStream << setprecision(precision)
                         << *multiply(dx / dt,
                                      *computeVelocityNorm(
                                          lattice, lattice.getBoundingBox()));
  // velocity components
  const string velocityFileName =
      outputDir + createFileName("velocity", iter, maxdigit) + ".dat";
  plb_ofstream velocityFileStream(velocityFileName.c_str());
  velocityFileStream << setprecision(precision)
                     << *multiply(dx / dt,
                                  *computeVelocity(lattice,
                                                   lattice.getBoundingBox()));
}

int main(int argc, char *argv[]) {
  plbInit(&argc, &argv);

  pcout << "Read Configuration..." << endl;
  XMLreader xmlFile("guo_case.xml");

  string outputDir;
  plint resolution;
  T lx;
  T ly;
  T re;
  T omega;
  T logT;
  T outputT;
  T maxT;
  T threshold;

  try {
    xmlFile["simuParam"]["re"].read(re);
    xmlFile["simuParam"]["omega"].read(omega);
    xmlFile["simuParam"]["resolution"].read(resolution);
    xmlFile["simuParam"]["threshold"].read(threshold);
    xmlFile["simuParam"]["lx"].read(lx);
    xmlFile["simuParam"]["ly"].read(ly);
    xmlFile["io"]["output"].read(outputDir);
    xmlFile["io"]["outputT"].read(outputT);
    xmlFile["io"]["logT"].read(logT);
    xmlFile["io"]["maxT"].read(maxT);
  } catch (const PlbIOException &e) {
    std::cerr << e.what() << endl;
    return -1;
  }
  global::directories().setOutputDir(outputDir);

  const T tau = 1.0 / omega;
  const T visco = (tau - 0.5) / 3.0;
  const plint diameter = 30;
  const T uo = re * visco / diameter;
  plint nx = lx * resolution + 1;
  plint ny = ly * resolution + 1;
  plint cx = nx / 14 * 2.5; // cx=75
  plint cy = ny / 2;        // cy=105
  pcout << "cx: " << cx << " ,cy: " << cy << " ,diameter: " << diameter << endl;
  T pl_re = resolution * uo / visco;

  IncomprFlowParam<T> parameters(uo,         // uMax
                                 pl_re,      // Re
                                 resolution, // N
                                 lx,         // lx
                                 ly          // ly
  );

  const plint logStep = parameters.nStep(logT);
  const plint residualAnalysisStep = parameters.nStep(logT);
  const plint imSaveStep = parameters.nStep(outputT);
  const plint vtkSaveStep = parameters.nStep(outputT);
  const plint fieldSaveStep = parameters.nStep(outputT);

  writeLogFile(parameters, "Poiseuille flow");

  MultiBlockLattice2D<T, DESCRIPTOR> lattice(
      parameters.getNx(), parameters.getNy(),
      new BGKdynamics<T, DESCRIPTOR>(parameters.getOmega()));

  OnLatticeBoundaryCondition2D<T, DESCRIPTOR> *boundaryCondition =
      createLocalBoundaryCondition2D<T, DESCRIPTOR>();

  cylinderSetup(lattice, parameters, *boundaryCondition, cx, cy, diameter / 2);

  util::ResidualTracer2D<T> residualTracer(10, nx, ny, threshold);
  MultiScalarField2D<T> previousVelNorm(nx, ny);
  MultiScalarField2D<T> currentVelNorm(nx, ny);
  Box2D domain = lattice.getBoundingBox();

  global::timer("mainloop").start();
  // Main loop over time iterations.
  plint iT = 0;
  for (iT = 0; iT * parameters.getDeltaT() < maxT; ++iT) {

    // At this point, the state of the lattice corresponds to the
    //   discrete time iT. However, the stored averages
    //   (getStoredAverageEnergy and getStoredAverageDensity) correspond to
    //   the previous time iT-1.

    if (iT % residualAnalysisStep == residualAnalysisStep - 1) {
      computeVelocityNorm(lattice, previousVelNorm, domain);
    }

    if (iT % residualAnalysisStep == 0 && iT > 0) {
      computeVelocityNorm(lattice, currentVelNorm, domain);
      residualTracer.measure(currentVelNorm, previousVelNorm, domain, true);
    }

    if (residualTracer.hasConverged()) {
      pcout << "simulation is over" << endl;
      break;
    }

    // Lattice Boltzmann iteration step.
    lattice.collideAndStream();

    // At this point, the state of the lattice corresponds to the
    //   discrete time iT+1, and the stored averages are upgraded to time iT.
    if (iT % logStep == 0) {
      pcout << "step " << iT << "; t=" << iT * parameters.getDeltaT()
            << "; av energy =" << setprecision(10)
            << getStoredAverageEnergy<T>(lattice)
            << "; av rho =" << getStoredAverageDensity<T>(lattice) << endl;
    }

    if (iT % imSaveStep == 0) {
      pcout << "Saving Gif ..." << endl;
      writeGif(lattice, iT);
    }

    if (iT % vtkSaveStep == 0 && iT > 0) {
      pcout << "Saving VTK file ..." << endl;
      writeVTK(lattice, parameters, iT);
    }

    if (iT % fieldSaveStep == 0 && iT > 0) {
      pcout << "Saving field text file..." << endl;
      writeField(outputDir, lattice, parameters, iT);
    }
  }

  pcout << "write out the fields at the final step..." << endl;
  writeGif(lattice, iT);
  writeVTK(lattice, parameters, iT);
  writeField(outputDir, lattice, parameters, iT);

  T tEnd = global::timer("mainloop").stop();
  pcout << "number of processors: " << global::mpi().getSize() << endl;
  pcout << "total iteraction step: " << iT << endl;
  pcout << "total cpu clock time for main loop: " << tEnd << endl;

  delete boundaryCondition;
}
