#include "../../src/parameters/singlePhaseFlowParameters.h"
#include "basicDynamics/isoThermalDynamics.h"
#include "core/globalDefs.h"
#include "core/plbInit.h"
#include "core/runTimeDiagnostics.h"
#include "io/imageWriter.h"
#include "io/parallelIO.h"
#include "libraryInterfaces/TINYXML_xmlIO.h"
#include "multiBlock/multiBlockLattice2D.h"
#include "multiBlock/multiDataField2D.h"
#include "palabos2D.h"
#include "palabos2D.hh"
#include <memory>

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
  Param() = default;
  Param(string configXmlName) {
    XMLreader document(configXmlName);
    document["geometry"]["filename"].read(geometryFile);
    document["geometry"]["nx"].read(nx);
    document["geometry"]["ny"].read(ny);
    document["io"]["output"].read(outputDir);
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
    pcout << "geometry image is output: " << param.outputDir << seperator
          << imageFileName << ".gif" << endl;
  }
  return geometry;
}

int main(int argc, char **argv) {
  plbInit(&argc, &argv);
  string configXml;
  try {
    global::argv(1).read(configXml);

  } catch (const PlbIOException &ex) {
    pcerr << "Wrong parameters; the syntax is" << (string)global::argv(0)
          << " config.xml" << endl;
  }
  Param param;
  try {
    param = Param(configXml);
  } catch (const PlbIOException &ex) {
    pcerr << "Invaid configuration xml, with error message: " << ex.what()
          << endl;
  }
  pcout << "configuration xml file: " << configXml << endl;
  pcout << "geometry file: " << param.geometryFile << endl;
  pcout << "output directory: " << param.outputDir << endl;
  global::directories().setOutputDir(param.outputDir);
  param.flowParam.writeLogFile("Berea Stone 2D flow");

  std::unique_ptr<MultiScalarField2D<int>> geometry = readGeomtry(param);
  MultiBlockLattice2D<T, DESCRIPTOR> lattices(
      param.nx, param.ny,
      new BGKdynamics<T, DESCRIPTOR>(param.flowParam.getOmega()));

  return 0;
}