/** \file
 * Flow around a 2D cylinder inside a channel, with the creation of a von
 * Karman vortex street. This example makes use of bounce-back nodes to
 * describe the shape of the cylinder. The outlet is modeled through a
 * Neumann (zero velocity-gradient) condition.
 */

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

/// Velocity on the parabolic Poiseuille profile
T poiseuilleVelocity(plint iY, IncomprFlowParam<T> const &parameters) {
  T y = (T)iY / parameters.getResolution();
  return 4. * parameters.getLatticeU() * (y - y * y);
}

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

/// A functional, used to initialize the velocity for the boundary conditions
template <typename T> class PoiseuilleVelocity {
public:
  PoiseuilleVelocity(IncomprFlowParam<T> parameters_)
      : parameters(parameters_) {}
  void operator()(plint iX, plint iY, Array<T, 2> &u) const {
    u[0] = poiseuilleVelocity(iY, parameters);
    u[1] = T();
  }

private:
  IncomprFlowParam<T> parameters;
};

/// A functional, used to initialize a pressure boundary to constant density
template <typename T> class ConstantDensity {
public:
  ConstantDensity(T density_) : density(density_) {}
  T operator()(plint iX, plint iY) const { return density; }

private:
  T density;
};

/// A functional, used to create an initial condition for the density and
/// velocity
template <typename T> class PoiseuilleVelocityAndDensity {
public:
  PoiseuilleVelocityAndDensity(IncomprFlowParam<T> parameters_)
      : parameters(parameters_) {}
  void operator()(plint iX, plint iY, T &rho, Array<T, 2> &u) const {
    rho = poiseuilleDensity(iX, parameters);
    u[0] = poiseuilleVelocity(iY, parameters);
    u[1] = T();
  }

private:
  IncomprFlowParam<T> parameters;
};

/// A functional, used to instantiate bounce-back nodes at the locations of the
/// cylinder
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
    OnLatticeBoundaryCondition2D<T, DESCRIPTOR> &boundaryCondition) {
  const plint nx = parameters.getNx();
  const plint ny = parameters.getNy();
  Box2D outlet(nx - 1, nx - 1, 1, ny - 2);

  // Create Velocity boundary conditions everywhere
  boundaryCondition.setVelocityConditionOnBlockBoundaries(
      lattice, Box2D(0, 0, 1, ny - 2));
  boundaryCondition.setVelocityConditionOnBlockBoundaries(
      lattice, Box2D(0, nx - 1, 0, 0));
  boundaryCondition.setVelocityConditionOnBlockBoundaries(
      lattice, Box2D(0, nx - 1, ny - 1, ny - 1));
  // .. except on right boundary, where we prefer an outflow condition
  //    (zero velocity-gradient).
  boundaryCondition.setVelocityConditionOnBlockBoundaries(
      lattice, Box2D(nx - 1, nx - 1, 1, ny - 2), boundary::outflow);

  setBoundaryVelocity(lattice, lattice.getBoundingBox(),
                      PoiseuilleVelocity<T>(parameters));
  setBoundaryDensity(lattice, outlet, ConstantDensity<T>(1.));
  initializeAtEquilibrium(lattice, lattice.getBoundingBox(),
                          PoiseuilleVelocityAndDensity<T>(parameters));

  plint cx = nx / 4;
  plint cy = ny / 2 + 2; // cy is slightly offset to avoid full symmetry,
                         //   and to get a Von Karman Vortex street.
  plint radius = cy / 4;
  defineDynamics(lattice, lattice.getBoundingBox(),
                 new CylinderShapeDomain2D<T>(cx, cy, radius),
                 new plb::BounceBack<T, DESCRIPTOR>);

  lattice.initialize();
}

void writeGif(MultiBlockLattice2D<T, DESCRIPTOR> &lattice, plint iter) {
  ImageWriter<T> imageWriter("leeloo");
  imageWriter.writeScaledGif(createFileName("u", iter, 6),
                             *computeVelocityNorm(lattice));
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
                IncomprFlowParam<T> const &parameters, plint iter) {
  T dx = parameters.getDeltaX();
  T dt = parameters.getDeltaT();
  int precision = 6;
  // density
  const string densityFileName =
      outputDir + createFileName("density", iter, 6) + ".dat";
  plb_ofstream densityFileStream(densityFileName.c_str());
  MultiScalarField2D<T> densityField =
      *computeDensity(lattice, lattice.getBoundingBox());
  densityFileStream << setprecision(precision) << densityField;
  // velocity norm
  const string velocityNormFileName =
      outputDir + createFileName("velocityNorm", iter, 6) + ".dat";
  plb_ofstream velocityNormFileStream(velocityNormFileName.c_str());
  MultiScalarField2D<T> velocityNormField = *multiply(
      dx / dt, *computeVelocityNorm(lattice, lattice.getBoundingBox()));
  velocityNormFileStream << setprecision(precision) << velocityNormField;
  // velocity components
  const string velocityFileName =
      outputDir + createFileName("velocity", iter, 6) + ".dat";
  plb_ofstream velocityFileStream(velocityFileName.c_str());
  MultiTensorField2D<T, 2> velocityField =
      *multiply(dx / dt, *computeVelocity(lattice, lattice.getBoundingBox()));
  velocityFileStream << setprecision(precision) << velocityField;
}

int main(int argc, char *argv[]) {
  plbInit(&argc, &argv);

  const string outputDir = "./tmp/";
  global::directories().setOutputDir(outputDir);

  IncomprFlowParam<T> parameters((T)1e-2, // uMax
                                 (T)600., // Re
                                 100,     // N
                                 6.,      // lx
                                 1.       // ly
  );
  const T logT = (T)0.02;
  const T imSave = (T)0.06;
  const T vtkSave = (T)1.;
  const T filedSave = T(1.);
  const T maxT = (T)3.1;

  writeLogFile(parameters, "Poiseuille flow");

  MultiBlockLattice2D<T, DESCRIPTOR> lattice(
      parameters.getNx(), parameters.getNy(),
      new BGKdynamics<T, DESCRIPTOR>(parameters.getOmega()));


/*
	for (int i = (-20+parameters.getNx()/2); i < (20+parameters.getNx()/2); ++i)
	{
		for (int j = (-20+parameters.getNy()/2); j < (20+parameters.getNy()/2); ++j)
		{
			if (pow(i-50,2)+pow(j-50,2)<=400)
			    lattice.get(i,j).getDynamics().setOmega(1);
			else
			    lattice.get(i,j).getDynamics().setOmega(0.5);
		}
	}*/

  OnLatticeBoundaryCondition2D<T, DESCRIPTOR> *boundaryCondition =
      createLocalBoundaryCondition2D<T, DESCRIPTOR>();

  cylinderSetup(lattice, parameters, *boundaryCondition);

  // Main loop over time iterations.
  for (plint iT = 0; iT * parameters.getDeltaT() < maxT; ++iT) {
    // At this point, the state of the lattice corresponds to the
    //   discrete time iT. However, the stored averages (getStoredAverageEnergy
    //   and getStoredAverageDensity) correspond to the previous time iT-1.

    // if (iT % parameters.nStep(imSave) == 0) {
    //   pcout << "Saving Gif ..." << endl;
    //   writeGif(lattice, iT);
    // }

    if (iT % parameters.nStep(vtkSave) == 0 && iT > 0) {
      pcout << "Saving VTK file ..." << endl;
      writeVTK(lattice, parameters, iT);
    }

    if (iT % parameters.nStep(filedSave) == 0 && iT >= 0) {
      pcout << "Saving filed text file..." << endl;
      writeField(outputDir, lattice, parameters, iT);
    }

    if (iT % parameters.nStep(logT) == 0) {
      pcout << "step " << iT << "; t=" << iT * parameters.getDeltaT();
    }

    // Lattice Boltzmann iteration step.
    lattice.collideAndStream();

    // At this point, the state of the lattice corresponds to the
    //   discrete time iT+1, and the stored averages are upgraded to time iT.
    if (iT % parameters.nStep(logT) == 0) {
      pcout << "; av energy =" << setprecision(10)
            << getStoredAverageEnergy<T>(lattice)
            << "; av rho =" << getStoredAverageDensity<T>(lattice) << endl;
    }
  }

  delete boundaryCondition;
}
