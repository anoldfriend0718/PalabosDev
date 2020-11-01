#pragma once

#include "customizedutil/utilities.h"
#include "io/imageWriter.h"
#include "io/vtkDataOutput.h"
#include "multiBlock/multiBlockLattice2D.h"
#include "multiBlock/multiDataField2D.h"
#include "parameters/LBMModelParser2D.h"
#include "parameters/singlePhaseFlowParameters.h"

namespace plb {
namespace util {

namespace sp {

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField2D<int>>
readGeomtry(const SinglePhaseFlowCaseParam<T, Descriptor> &param,
            bool isWriteGeometryImage = true) {
  std::unique_ptr<MultiScalarField2D<int>> geometry(
      new MultiScalarField2D<int>(param.nx, param.ny));
  plb_ifstream geometryFileStream(param.geometryFile.c_str());
  geometryFileStream >> *geometry;
  if (isWriteGeometryImage) {
    ImageWriter<int> geometryImageWriter("leeloo");
    const std::string seperator = "/";
    std::string imageFileName =
        util::getFileName(param.geometryFile, seperator);
    geometryImageWriter.writeScaledGif(imageFileName, *geometry);
    PLOG(plog::info) << "geometry image is output: " << param.outputDir
                     << seperator << imageFileName << ".gif";
  }
  return geometry;
}

template <typename T, template <typename U> class Descriptor>
void writeGif(MultiBlockLattice2D<T, Descriptor> &lattice, plint iter,
              plint maxdigit = 8) {
  ImageWriter<T> imageWriter("leeloo");
  imageWriter.writeScaledGif(createFileName("u", iter, maxdigit),
                             *computeVelocityNorm(lattice));
  imageWriter.writeScaledGif(createFileName("rho", iter, maxdigit),
                             *computeDensity(lattice));
}

template <typename T, template <typename U> class Descriptor>
void writeVTK(MultiBlockLattice2D<T, Descriptor> &lattice,
              SinglePhaseFlowParam<T> const &parameters, plint iter,
              plint maxdigit = 8) {
  T dx = parameters.getDeltaX();
  T dt = parameters.getDeltaT();
  VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, maxdigit), dx);
  vtkOut.template writeData<float>(*computeVelocityNorm(lattice),
                                   "velocityNorm", dx / dt);
  vtkOut.template writeData<2, float>(*computeVelocity(lattice), "velocity",
                                      dx / dt);
}

template <typename T, template <typename U> class Descriptor>
void writeField(const std::string outputDir,
                MultiBlockLattice2D<T, Descriptor> &lattice,
                SinglePhaseFlowParam<T> const &parameters, plint iter,
                plint maxdigit = 8) {
  T dx = parameters.getDeltaX();
  T dt = parameters.getDeltaT();
  int precision = 6;
  // density
  const std::string densityFileName =
      outputDir + createFileName("density", iter, maxdigit) + ".dat";
  plb_ofstream densityFileStream(densityFileName.c_str());
  densityFileStream << std::setprecision(precision)
                    << *computeDensity(lattice, lattice.getBoundingBox());
  // velocity norm
  const std::string velocityNormFileName =
      outputDir + createFileName("velocityNorm", iter, maxdigit) + ".dat";
  plb_ofstream velocityNormFileStream(velocityNormFileName.c_str());
  velocityNormFileStream << std::setprecision(precision)
                         << *multiply(dx / dt,
                                      *computeVelocityNorm(
                                          lattice, lattice.getBoundingBox()));
  // velocity components
  const std::string velocityFileName =
      outputDir + createFileName("velocity", iter, maxdigit) + ".dat";
  plb_ofstream velocityFileStream(velocityFileName.c_str());
  velocityFileStream << std::setprecision(precision)
                     << *multiply(dx / dt,
                                  *computeVelocity(lattice,
                                                   lattice.getBoundingBox()));
}

template <typename T, template <typename U> class Descriptor>
T computePermeability(MultiBlockLattice2D<T, Descriptor> &lattice,
                      SinglePhaseFlowParam<T> const &flowParam) {
  Box2D domain = lattice.getBoundingBox();
  pluint nx = flowParam.getNx();
  pluint ny = flowParam.getNy();
  Box2D inlet(0, 0, 1, ny - 2);
  Box2D outlet(nx - 1, nx - 1, 1, ny - 2);
  T densityInlet = computeAverageDensity(lattice, inlet);
  T densityOutlet = computeAverageDensity(lattice, outlet);
  T deltaP = (densityInlet - densityOutlet) * Descriptor<T>::cs2;
  PLOG(plog::debug) << "inlet average density: " << densityInlet
                    << "; outlet average density: " << densityOutlet
                    << "; delta pressure: " << deltaP;
  std::unique_ptr<MultiScalarField2D<T>> velU =
      computeVelocityComponent(lattice, domain, 0);
  T meanUInlet = computeAverage(*velU, inlet);
  T meanUOutlet = computeAverage(*velU, outlet);
  T meanUDomain = computeAverage(*velU, domain);
  PLOG(plog::debug) << "inlet average U: " << meanUInlet
                    << "; outlet average U: " << meanUOutlet
                    << "; domain average U: " << meanUDomain;
  T permeability =
      flowParam.getLatticeNu() * meanUDomain / (deltaP / (T)(nx - 1));
  PLOG(plog::info) << "permeability: " << permeability;
  return permeability;
}

/// Velocity on the parabolic Poiseuille profile
template <typename T>
T poiseuilleVelocityProfile(plint iY,
                            SinglePhaseFlowParam<T> const &parameters) {
  T y = (T)iY / parameters.getNy();
  return 4. * parameters.getLatticeU() * (y - y * y);
}

/// A functional, used to initialize the velocity for the boundary conditions
template <typename T> class PoiseuilleVelocity {
public:
  PoiseuilleVelocity(SinglePhaseFlowParam<T> parameters_)
      : parameters(parameters_) {}
  void operator()(plint iX, plint iY, Array<T, 2> &u) const {
    u[0] = poiseuilleVelocityProfile(iY, parameters);
    u[1] = T();
  }

private:
  SinglePhaseFlowParam<T> parameters;
};

} // namespace sp

} // namespace util
} // namespace plb