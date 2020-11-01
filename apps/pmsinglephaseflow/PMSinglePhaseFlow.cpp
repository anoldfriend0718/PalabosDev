
#include "basicDynamics/isoThermalDynamics.h"
#include "boundaryCondition/boundaryCondition.h"
#include "core/globalDefs.h"
#include "core/plbInit.h"
#include "core/plbProfiler.h"
#include "core/runTimeDiagnostics.h"
#include "customizedutil/plblogger.h"
#include "customizedutil/residualTracer.h"
#include "customizedutil/residualTracer.hh"
#include "customizedutil/singlePhaseFlowUtities.h"
#include "customizedutil/utilities.h"
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
#include "parameters/LBMModelParser2D.h"
#include "parameters/LBMModelParser2D.hh"
#include "parameters/singlePhaseFlowParameters.h"
#include "plog/Log.h"
#include "plog/Severity.h"
#include <cstdlib>
#include <memory>
#include <plog/Formatters/TxtFormatter.h>

using namespace plb;
using namespace plb::util::sp;
using namespace std;

#define DESCRIPTOR descriptors::D2Q9Descriptor
typedef double T;

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
                      plb::util::sp::PoiseuilleVelocity<T>(flowParam));

  // Where "geometry" evaluates to 1, use bounce-back.
  defineDynamics(lattice, *geometry, new BounceBack<T, DESCRIPTOR>(), 1);
  // Where "geometry" evaluates to 2, use no-dynamics (which does nothing).
  defineDynamics(lattice, *geometry, new NoDynamics<T, DESCRIPTOR>(), 2);

  initializeAtEquilibrium(lattice, lattice.getBoundingBox(),
                          ZeroVelocityAndConstDensity<T>(flowParam));
  lattice.initialize();
  PLOG(plog::info) << "simulation is set up";
}

std::string
getCheckPointFile(const SinglePhaseFlowCaseParam<T, DESCRIPTOR> &param,
                  pluint timeStep, pluint maxDigit = 8) {
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
  SinglePhaseFlowCaseParam<T, DESCRIPTOR> param;
  try {
    param = SinglePhaseFlowCaseParam<T, DESCRIPTOR>(configXml);
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

  std::unique_ptr<Dynamics<T, DESCRIPTOR>> dynamics =
      param.lbmPara.getDynamics();
  PLOG(plog::info) << "collision dynamics: "
                   << util::getDynamicsName<T, DESCRIPTOR>(dynamics->getId());
  // lattice help release the dnamics pointer, therefore we need to clone
  // dynamics object to avoid duplicated deleting same pointer
  MultiBlockLattice2D<T, DESCRIPTOR> lattice(param.nx, param.ny,
                                             dynamics->clone());
  boundarySetAndInit(lattice, param.flowParam, geometry);

  global::timer("mainloop").start();
  global::profiler().turnOn();

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
  for (; iT * deltaT <= param.maxT; iT++) {
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
  global::profiler().writeReport();
  return EXIT_SUCCESS;
}