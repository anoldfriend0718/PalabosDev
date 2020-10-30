#include "../../src/parameters/singlePhaseFlowParameters.h"
#include "basicDynamics/isoThermalDynamics.h"
#include "boundaryCondition/boundaryCondition.h"
#include "core/globalDefs.h"
#include "core/plbInit.h"
#include "core/runTimeDiagnostics.h"
#include "customizedutil/plblogger.h"
#include "customizedutil/residualTracer.h"
#include "customizedutil/residualTracer.hh"
#include "dataProcessors/dataAnalysisWrapper2D.h"
#include "dataProcessors/dataAnalysisWrapper3D.h"
#include "dataProcessors/dataInitializerWrapper2D.h"
#include "io/imageWriter.h"
#include "io/parallelIO.h"
#include "io/serializerIO_2D.h"
#include "libraryInterfaces/TINYXML_xmlIO.h"
#include "multiBlock/multiBlockLattice2D.h"
#include "multiBlock/multiDataField2D.h"
#include "palabos2D.h"
#include "palabos2D.hh"
#include "plog/Log.h"
#include "plog/Severity.h"
#include <cstdlib>
#include <memory>
#include <plog/Formatters/TxtFormatter.h>

using namespace plb;
using namespace std;

#define DESCRIPTOR descriptors::D2Q9Descriptor
typedef double T;

struct Param {
  string geometryFile;
  string outputDir;
  pluint nx;
  pluint ny;
  SinglePhaseFlowParam<T> flowParam;
  T threshold;
  pluint logStep;
  pluint imSaveStep;
  pluint vtkSaveStep;
  pluint fieldDataSaveStep;
  pluint residualAnalysisStep;
  T maxT;
  std::string logLevel;
  bool ifImSave;
  bool ifVtkSave;
  bool ifFiledDataSave;
  bool ifToggleInternalStatistics;
  bool ifCheckPoint;
  string checkPointFilePrefix;
  pluint checkPointStep;
  bool ifLoadCheckPoint;
  pluint initStep;
  Param() = default;
  Param(string configXmlName) {
    XMLreader document(configXmlName);
    document["geometry"]["filename"].read(geometryFile);
    document["geometry"]["nx"].read(nx);
    document["geometry"]["ny"].read(ny);
    T physicalU;
    pluint latticeCharacteristicLength;
    T resolution;
    T re;
    T latticeU;
    pluint latticeLx;
    pluint latticeLy;
    document["simuParam"]["physicalU"].read(physicalU);
    document["simuParam"]["latticeU"].read(latticeU);
    document["simuParam"]["re"].read(re);
    document["simuParam"]["latticeLx"].read(latticeLx);
    document["simuParam"]["latticeLy"].read(latticeLy);
    document["simuParam"]["resolution"].read(resolution);
    document["simuParam"]["latticeCharacteristicLength"].read(
        latticeCharacteristicLength);
    flowParam = SinglePhaseFlowParam<T>(physicalU, latticeU, re,
                                        latticeCharacteristicLength, resolution,
                                        latticeLx, latticeLy);
    document["simuParam"]["resolution"].read(threshold);
    document["io"]["output"].read(outputDir);
    document["io"]["logStep"].read(logStep);
    document["io"]["imSaveStep"].read(imSaveStep);
    document["io"]["vtkSaveStep"].read(vtkSaveStep);
    document["io"]["fieldDataSaveStep"].read(fieldDataSaveStep);
    document["io"]["residualAnalysisStep"].read(residualAnalysisStep);
    document["io"]["maxT"].read(maxT);
    document["io"]["logLevel"].read(logLevel);
    document["io"]["ifImSave"].read(ifImSave);
    document["io"]["ifVtkSave"].read(ifVtkSave);
    document["io"]["ifFiledDataSave"].read(ifFiledDataSave);
    document["io"]["ifToggleInternalStatistics"].read(
        ifToggleInternalStatistics);
    document["io"]["ifCheckPoint"].read(ifCheckPoint);
    document["io"]["checkPointFilePrefix"].read(checkPointFilePrefix);
    document["io"]["checkPointStep"].read(checkPointStep);
    document["init"]["ifLoadCheckPoint"].read(ifLoadCheckPoint);
    document["init"]["initStep"].read(initStep);
  }
};

std::string GetFileName(const std::string path, const string seperator) {
  string::size_type iPos = path.find_last_of(seperator) + 1;
  string filename = path.substr(iPos, path.length() - iPos);
  string name = filename.substr(0, filename.rfind("."));
  return name;
}

unique_ptr<MultiScalarField2D<int>>
readGeomtry(const Param &param, bool isWriteGeometryImage = true) {
  unique_ptr<MultiScalarField2D<int>> geometry(
      new MultiScalarField2D<int>(param.nx, param.ny));
  plb_ifstream geometryFileStream(param.geometryFile.c_str());
  geometryFileStream >> *geometry;
  if (isWriteGeometryImage) {
    ImageWriter<int> geometryImageWriter("leeloo");
    const string seperator = "/";
    string imageFileName = GetFileName(param.geometryFile, seperator);
    geometryImageWriter.writeScaledGif(imageFileName, *geometry);
    PLOG(plog::info) << "geometry image is output: " << param.outputDir
                     << seperator << imageFileName << ".gif" << endl;
  }
  return geometry;
}

/// Velocity on the parabolic Poiseuille profile
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

template <typename T> class ZeroVelocityAndConstDensity {
public:
  ZeroVelocityAndConstDensity(SinglePhaseFlowParam<T> parameters_)
      : parameters(parameters_) {}
  void operator()(plint iX, plint iY, T &rho, Array<T, 2> &u) const {
    rho = 1.0;
    u[0] = 0.0;
    u[1] = T();
  }

private:
  SinglePhaseFlowParam<T> parameters;
};

void boundarySetAndInit(MultiBlockLattice2D<T, DESCRIPTOR> &lattice,
                        const SinglePhaseFlowParam<T> flowParam,
                        unique_ptr<MultiScalarField2D<int>> &geometry) {
  // top and bottom are set as periodic B.C.
  lattice.periodicity().toggle(1, true);
  // inlet and outlet are set as Regularized Dirichlet Velocity B.C.
  unique_ptr<OnLatticeBoundaryCondition2D<T, DESCRIPTOR>> boundaryCondition(
      createLocalBoundaryCondition2D<T, DESCRIPTOR>());
  pluint nx = flowParam.getNx();
  pluint ny = flowParam.getNy();
  Box2D inlet(0, 0, 1, ny - 2);
  Box2D outlet(nx - 1, nx - 1, 1, ny - 2);
  boundaryCondition->setVelocityConditionOnBlockBoundaries(lattice, inlet);
  boundaryCondition->setVelocityConditionOnBlockBoundaries(lattice, outlet,
                                                           boundary::outflow);
  setBoundaryVelocity(lattice, lattice.getBoundingBox(),
                      PoiseuilleVelocity<T>(flowParam));

  // Where "geometry" evaluates to 1, use bounce-back.
  defineDynamics(lattice, *geometry, new BounceBack<T, DESCRIPTOR>(), 1);
  // Where "geometry" evaluates to 2, use no-dynamics (which does nothing).
  defineDynamics(lattice, *geometry, new NoDynamics<T, DESCRIPTOR>(), 2);

  initializeAtEquilibrium(lattice, lattice.getBoundingBox(),
                          ZeroVelocityAndConstDensity<T>(flowParam));
  lattice.initialize();
  PLOG(plog::info) << "simulation is set up";
}

void writeGif(MultiBlockLattice2D<T, DESCRIPTOR> &lattice, plint iter,
              plint maxdigit = 8) {
  ImageWriter<T> imageWriter("leeloo");
  imageWriter.writeScaledGif(createFileName("u", iter, maxdigit),
                             *computeVelocityNorm(lattice));
  imageWriter.writeScaledGif(createFileName("rho", iter, maxdigit),
                             *computeDensity(lattice));
}

void writeVTK(MultiBlockLattice2D<T, DESCRIPTOR> &lattice,
              SinglePhaseFlowParam<T> const &parameters, plint iter,
              plint maxdigit = 8) {
  T dx = parameters.getDeltaX();
  T dt = parameters.getDeltaT();
  VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, maxdigit), dx);
  vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm",
                          dx / dt);
  vtkOut.writeData<2, float>(*computeVelocity(lattice), "velocity", dx / dt);
}

void writeField(const string outputDir,
                MultiBlockLattice2D<T, DESCRIPTOR> &lattice,
                SinglePhaseFlowParam<T> const &parameters, plint iter,
                plint maxdigit = 8) {
  T dx = parameters.getDeltaX();
  T dt = parameters.getDeltaT();
  int precision = 6;
  // density
  const string densityFileName =
      outputDir + createFileName("density", iter, maxdigit) + ".dat";
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

T computePermeability(MultiBlockLattice2D<T, DESCRIPTOR> &lattice,
                      SinglePhaseFlowParam<T> const &flowParam) {
  Box2D domain = lattice.getBoundingBox();
  pluint nx = flowParam.getNx();
  pluint ny = flowParam.getNy();
  Box2D inlet(0, 0, 1, ny - 2);
  Box2D outlet(nx - 1, nx - 1, 1, ny - 2);
  T densityInlet = computeAverageDensity(lattice, inlet);
  T densityOutlet = computeAverageDensity(lattice, outlet);
  T deltaP = (densityInlet - densityOutlet) * DESCRIPTOR<T>::cs2;
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

std::string getCheckPointFile(const Param &param, pluint timeStep,
                              pluint maxDigit = 8) {
  std::string outputDir = param.outputDir;
  std::string filePath =
      outputDir + "/" +
      createFileName(param.checkPointFilePrefix, timeStep, maxDigit);
  return filePath;
}

int main(int argc, char **argv) {
  plbInit(&argc, &argv);

  string configXml;
  try {
    global::argv(1).read(configXml);

  } catch (const PlbIOException &ex) {
    pcerr << "Wrong parameters; the syntax is" << (string)global::argv(0)
          << " config.xml" << endl;
    return EXIT_FAILURE;
  }
  Param param;
  try {
    param = Param(configXml);
  } catch (const PlbIOException &ex) {
    pcerr << "Invaid configuration xml, with error message: " << ex.what()
          << endl;
    return EXIT_FAILURE;
  }
  // initialize the logger
  static plb::util::PlbLoggerAppender<plog::TxtFormatter> plbAppender;
  plog::init(plb::util::getSeverity(param.logLevel), &plbAppender);

  PLOG(plog::info) << "configuration xml file: " << configXml;
  PLOG(plog::info) << "geometry file: " << param.geometryFile;
  PLOG(plog::info) << "output directory: " << param.outputDir;
  global::directories().setOutputDir(param.outputDir);
  param.flowParam.writeLogFile("Berea Stone 2D flow");

  std::unique_ptr<MultiScalarField2D<int>> geometry = readGeomtry(param);
  MultiBlockLattice2D<T, DESCRIPTOR> lattice(
      param.nx, param.ny,
      new BGKdynamics<T, DESCRIPTOR>(param.flowParam.getOmega()));
  boundarySetAndInit(lattice, param.flowParam, geometry);

  global::timer("mainloop").start();
  pluint nx = param.flowParam.getNx();
  pluint ny = param.flowParam.getNy();
  util::ResidualTracer2D<T> residualTracer(1, nx, ny, param.threshold);
  MultiScalarField2D<T> previousVelNorm(nx, ny);
  MultiScalarField2D<T> currentVelNorm(nx, ny);
  Box2D domain = lattice.getBoundingBox();
  T deltaT = param.flowParam.getDeltaT();
  pluint residualAnalysisStep = param.residualAnalysisStep;

  pluint initStep = 0;
  if (param.ifLoadCheckPoint) {
    PLOG(plog::info) << "loading checking point data from time step: "
                     << param.initStep;
    loadBinaryBlock(lattice, getCheckPointFile(param, param.initStep));
    initStep = param.initStep;
  }
  pluint iT = initStep;

  PLOG(plog::info) << "starting main loop...";
  for (; iT * deltaT < param.maxT; iT++) {
    lattice.toggleInternalStatistics(param.ifToggleInternalStatistics);
    if (param.ifToggleInternalStatistics && iT % param.logStep == 0) {
      PLOG(plog::info) << "step " << iT
                       << "; t=" << iT * param.flowParam.getDeltaT()
                       << "; av energy =" << setprecision(10)
                       << getStoredAverageEnergy<T>(lattice)
                       << "; av rho =" << getStoredAverageDensity<T>(lattice);
    }

    if (param.ifCheckPoint && iT % param.checkPointStep == 0 && iT > initStep) {
      PLOG(plog::info) << "step " << iT
                       << "; t=" << iT * param.flowParam.getDeltaT()
                       << "; checking point ...";
      saveBinaryBlock(lattice, getCheckPointFile(param, iT));
    }

    if (param.ifImSave && iT % param.imSaveStep == 0) {
      PLOG(plog::debug) << "step " << iT
                        << "; t=" << iT * param.flowParam.getDeltaT()
                        << "; Saving Gif ...";
      writeGif(lattice, iT);
    }

    if (param.ifVtkSave && iT % param.vtkSaveStep == 0 && iT > initStep) {
      PLOG(plog::debug) << "step " << iT
                        << "; t=" << iT * param.flowParam.getDeltaT()
                        << "; Saving VTK file ...";
      writeVTK(lattice, param.flowParam, iT);
    }

    if (param.ifFiledDataSave && iT % param.fieldDataSaveStep == 0 &&
        iT > initStep) {
      PLOG(plog::debug) << "step " << iT
                        << "; t=" << iT * param.flowParam.getDeltaT()
                        << "; Saving field text file...";
      writeField(param.outputDir, lattice, param.flowParam, iT);
    }

    if (iT % residualAnalysisStep == residualAnalysisStep - 1 &&
        iT > initStep) {
      computeVelocityNorm(lattice, previousVelNorm, domain);
    }

    if (iT % residualAnalysisStep == 0 && iT > initStep) {
      computeVelocityNorm(lattice, currentVelNorm, domain);
      PLOG(plog::info) << "step " << iT
                       << "; t=" << iT * param.flowParam.getDeltaT()
                       << "; measure residual error...";
      residualTracer.measure(currentVelNorm, previousVelNorm, domain, true);
    }

    if (residualTracer.hasConverged()) {
      break;
    }

    // Lattice Boltzmann iteration step.
    lattice.collideAndStream();

    // At this point, the state of the lattice corresponds to the
    //   discrete time iT+1, and the stored averages are upgraded to time iT.
  }

  PLOG(plog::info) << "step " << iT
                   << "; t=" << iT * param.flowParam.getDeltaT()
                   << "; write out the fields at the final step...";
  writeGif(lattice, iT);
  writeVTK(lattice, param.flowParam, iT);
  writeField(param.outputDir, lattice, param.flowParam, iT);

  computePermeability(lattice, param.flowParam);

  T tEnd = global::timer("mainloop").stop();
  PLOG(plog::info) << "number of processors: " << global::mpi().getSize();
  PLOG(plog::info) << "total iteraction step: " << iT;
  PLOG(plog::info) << "total computational time for main loop: " << tEnd;
  return EXIT_SUCCESS;
}