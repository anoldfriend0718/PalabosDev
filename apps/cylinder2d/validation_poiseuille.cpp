#include "palabos2D.h"
#include "palabos2D.hh"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "Guo_Dynamics_with_Porosity/Porosity_BGKdynamics.h"
#include "Guo_Dynamics_with_Porosity/Porosity_BGKdynamics.hh"
#include "Guo_Dynamics_with_Porosity/Porosity_dynamicsTemplates.h" 
#include "Guo_Dynamics_with_Porosity/Porosity_dynamicsTemplates2D.h"
#include "Guo_Dynamics_with_Porosity/Porosity_addDynamicParams.h"

#include "Guo_Dynamics_with_Porosity/Porosity_externalForceDynamics.h"
#include "Guo_Dynamics_with_Porosity/Porosity_externalForceDynamics.hh" 
#include "Guo_Dynamics_with_Porosity/Porosity_externalForceTemplates.h"
#include "Guo_Dynamics_with_Porosity/Porosity_externalForceTemplates2D.h"

using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
#define DESCRIPTOR D2Q9Descriptor

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


/// Velocity on the parabolic Poiseuille profile
T poiseuilleVelocity(plint iY, IncomprFlowParam<T> const& parameters) {
    T y = (T)iY / parameters.getResolution();
    return 4.*parameters.getLatticeU() * (y-y*y);
}

/// Linearly decreasing pressure profile
T poiseuillePressure(plint iX, IncomprFlowParam<T> const& parameters) {
    T Lx = parameters.getNx()-1;
    T Ly = parameters.getNy()-1;
    return 8.*parameters.getLatticeNu()*parameters.getLatticeU() / (Ly*Ly) * (Lx/(T)2-(T)iX);
}

/// Convert pressure to density according to ideal gas law
T poiseuilleDensity(plint iX, IncomprFlowParam<T> const& parameters) {
    return poiseuillePressure(iX,parameters)*DESCRIPTOR<T>::invCs2 + (T)1;
}

/// A functional, used to initialize the velocity for the boundary conditions
template<typename T>
class PoiseuilleVelocity {
public:
    PoiseuilleVelocity(IncomprFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    void operator()(plint iX, plint iY, Array<T,2>& u) const {
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = T();
    }
private:
    IncomprFlowParam<T> parameters;
};

/// A functional, used to initialize the density for the boundary conditions
template<typename T>
class PoiseuilleDensity {
public:
    PoiseuilleDensity(IncomprFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    T operator()(plint iX, plint iY) const {
        return poiseuilleDensity(iX,parameters);
    }
private:
    IncomprFlowParam<T> parameters;
};

/// A functional, used to create an initial condition for with zero velocity,
///   and linearly decreasing pressure.
template<typename T>
class PoiseuilleDensityAndZeroVelocity {
public:
    PoiseuilleDensityAndZeroVelocity(IncomprFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    void operator()(plint iX, plint iY, T& rho, Array<T,2>& u) const {
        rho = poiseuilleDensity(iX,parameters);
        u[0] = T();
        u[1] = T();
    }
private:
    IncomprFlowParam<T> parameters;
};

enum InletOutletT {pressure, velocity};

void channelSetup( MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
                   IncomprFlowParam<T> const& parameters,
                   OnLatticeBoundaryCondition2D<T,DESCRIPTOR>& boundaryCondition,
                   InletOutletT inletOutlet )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();

    // Note: The following approach illustrated here works only with boun-
    //   daries which are located on the outmost cells of the lattice. For
    //   boundaries inside the lattice, you need to use the version of
    //   "setVelocityConditionOnBlockBoundaries" which takes two Box2D
    //   arguments.

    // Velocity boundary condition on bottom wall. 
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
            lattice, Box2D(0, nx-1, 0, 0) );
    // Velocity boundary condition on top wall. 
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
            lattice, Box2D(0, nx-1, ny-1, ny-1) );

    // Pressure resp. velocity boundary condition on the inlet and outlet.
    if (inletOutlet == pressure) {
        // Note: pressure boundary conditions are currently implemented
        //   only for edges of the boundary, and not for corner nodes.
        boundaryCondition.setPressureConditionOnBlockBoundaries (
                lattice, Box2D(0,0, 1,ny-2) );
        boundaryCondition.setPressureConditionOnBlockBoundaries (
                lattice, Box2D(nx-1,nx-1, 1,ny-2) );
    }
    else {
        boundaryCondition.setVelocityConditionOnBlockBoundaries (
                lattice, Box2D(0,0, 1,ny-2) );
        boundaryCondition.setVelocityConditionOnBlockBoundaries (
                lattice, Box2D(nx-1,nx-1, 1,ny-2) );
    }

    // Define the value of the imposed density on all nodes which have previously been
    //   defined to be pressure boundary nodes.
    setBoundaryDensity (
            lattice, lattice.getBoundingBox(),
            PoiseuilleDensity<T>(parameters) );
    // Define the value of the imposed velocity on all nodes which have previously been
    //   defined to be velocity boundary nodes.
    setBoundaryVelocity (
            lattice, lattice.getBoundingBox(),
            PoiseuilleVelocity<T>(parameters) );
    // Initialize all cells at an equilibrium distribution, with a velocity and density
    //   value of the analytical Poiseuille solution.
    initializeAtEquilibrium (
           lattice, lattice.getBoundingBox(),
           PoiseuilleDensityAndZeroVelocity<T>(parameters) );

    lattice.initialize();
}

/// Produce a GIF snapshot of the velocity-norm.
void writeGif(MultiBlockLattice2D<T,DESCRIPTOR>& lattice, plint iter)
{
    const plint imSize = 600;

    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif(createFileName("u", iter, 6),
                               *computeVelocityNorm(lattice),
                               imSize, imSize );
}

/// Write the full velocity and the velocity-norm into a VTK file.
void writeVTK(MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
              IncomprFlowParam<T> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();
    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
    vtkOut.writeData<2,float>(*computeVelocity(lattice), "velocity", dx/dt);
}

T computeRMSerror ( MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
                    IncomprFlowParam<T> const& parameters )
{
    MultiTensorField2D<T,2> analyticalVelocity(lattice);
    setToFunction( analyticalVelocity, analyticalVelocity.getBoundingBox(),
                   PoiseuilleVelocity<T>(parameters) );
    MultiTensorField2D<T,2> numericalVelocity(lattice);
    computeVelocity(lattice, numericalVelocity, lattice.getBoundingBox());

           // Divide by lattice velocity to normalize the error
    return 1./parameters.getLatticeU() *
           // Compute RMS difference between analytical and numerical solution
               std::sqrt( computeAverage( *computeNormSqr(
                              *subtract(analyticalVelocity, numericalVelocity)
                          ) ) );
}

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);
    auto outputdir = "./tmp/";
    global::directories().setOutputDir(outputdir);

    IncomprFlowParam<T> parameters(
            (T) 2e-2,  // uMax
            (T) 5,    // Re:0.01_100,5
            60,        // N
            3.,        // lx
            1.         // ly 
    );
    T porosity = 0.1;//porosity can not = 0
    //T Da       = 1e-5;
    T KVC      = ((T)1/(T)3)* (parameters.getTau()-(T)0.5);//注意粘性系数与tau之间的关系
    T k0       = (T)1;
    //(150*Da*KVC*KVC*parameters.getRe()*parameters.getRe()*(1-porosity)*(1-porosity))/(parameters.getPhysicalU()*parameters.getPhysicalU()*(T)1.75*(T)1.75);//1;//.e-12;应该单位换算之后取格子单位
    //pcout << k0 << endl;

    const T logT     = (T)0.01;
    const T imSave   = (T)0.01;
    const T vtkSave  = (T)2.;
    const T maxT     = (T)15.1;
    // Change this variable to "pressure" if you prefer a pressure boundary
    //   condition with Poiseuille profile for the inlet and the outlet.
    const InletOutletT inletOutlet = velocity;

    writeLogFile(parameters, "Poiseuille flow");

    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
              parameters.getNx(), parameters.getNy(),
              //new Porosity_BGKdynamics<T,DESCRIPTOR>(parameters.getOmega(),0.9) );
              new Porosity_GuoExternalForceBGKdynamics<T,DESCRIPTOR>(parameters.getOmega(),porosity, KVC, k0) );
    
    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();

    channelSetup(lattice, parameters, *boundaryCondition, inletOutlet);

    // Main loop over time iterations.
    for (plint iT=0; iT*parameters.getDeltaT()<maxT; ++iT) {
       if (iT%parameters.nStep(imSave)==0) {
            pcout << "Saving Gif ..." << endl;
            writeGif(lattice, iT);
        }

        if (iT%parameters.nStep(vtkSave)==0 && iT>0) {
            pcout << "Saving VTK file ..." << endl;
            writeVTK(lattice, parameters, iT);
        }

        if (iT%parameters.nStep(imSave)==0 && iT>0) {
            pcout << "Saving field text file..." << endl;
            writeField(outputdir, lattice, parameters, iT);
        }

        if (iT%parameters.nStep(logT)==0) {
            pcout << "step " << iT
                  << "; t=" << iT*parameters.getDeltaT()
                  << "; RMS error=" << computeRMSerror(lattice, parameters);
            Array<T,2> uCenter;
            lattice.get(parameters.getNx()/2,parameters.getNy()/2).computeVelocity(uCenter);
            pcout << "; center velocity=" << uCenter[0]/parameters.getLatticeU() << endl;
        }

        // Lattice Boltzmann iteration step.
        lattice.collideAndStream();
    }

    delete boundaryCondition;
}
