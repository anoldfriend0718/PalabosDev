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

using namespace plb;
using namespace std;

#define DESCRIPTOR descriptors::D2Q9Descriptor
typedef double T;

struct Param {
  string geometryFile;
  string outputDir;
  pluint nx;
  pluint ny;
  Param() {}
  Param(string configXmlName) {
    XMLreader document(configXmlName);
    document["geometry"]["filename"].read(geometryFile);
    document["geometry"]["nx"].read(nx);
    document["geometry"]["ny"].read(ny);
    document["io"]["output"].read(outputDir);
  }
};

std::string GetFileName(const std::string path, const string seperator) {
  string::size_type iPos = path.find_last_of(seperator) + 1;
  string filename = path.substr(iPos, path.length() - iPos);
  string name = filename.substr(0, filename.rfind("."));
  return name;
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

  MultiScalarField2D<int> geometry(param.nx, param.ny);
  plb_ifstream geometryFileStream(param.geometryFile.c_str());
  geometryFileStream >> geometry;

  MultiBlockLattice2D<T, DESCRIPTOR> lattices(
      param.nx, param.ny, new BGKdynamics<T, DESCRIPTOR>(1.0));

  ImageWriter<int> geometryImageWriter("leeloo");
  geometryImageWriter.writeScaledGif(GetFileName(param.geometryFile, "/"),
                                     geometry);

  return 0;
}