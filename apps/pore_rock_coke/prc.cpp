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

/// A functional, used to initialize a pressure boundary to constant density
template<typename T>
class ConstantDensity {
public:
    ConstantDensity(T density_)
        : density(density_)
    { }
    T operator()(plint iX, plint iY) const {
        return density;
    }
private:
    T density;
};

/// A functional, used to create an initial condition for the density and velocity
template<typename T>
class PoiseuilleVelocityAndDensity {
public:
    PoiseuilleVelocityAndDensity(IncomprFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    void operator()(plint iX, plint iY, T& rho, Array<T,2>& u) const {
        rho = poiseuilleDensity(iX,parameters);
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = T();
    }
private:
    IncomprFlowParam<T> parameters;
};

void cylinderSetup( MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
                    IncomprFlowParam<T> const& parameters,
                    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>& boundaryCondition,
                    T omega,
                    T rockporosity,
                    T cokeporosity,
                    T KVC,
                    T k0 ,
                    MultiScalarField2D<int> boolMask)
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    Box2D outlet(nx-1,nx-1, 1, ny-2);

    // Create Velocity boundary conditions everywhere
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
            lattice, Box2D(0, 0, 1, ny-2) );
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
            lattice, Box2D(0, nx-1, 0, 0) );
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
            lattice, Box2D(0, nx-1, ny-1, ny-1) );
    // .. except on right boundary, where we prefer an outflow condition
    //    (zero velocity-gradient).
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
            lattice, Box2D(nx-1, nx-1, 1, ny-2), boundary::outflow );

    setBoundaryVelocity (
            lattice, lattice.getBoundingBox(),
            PoiseuilleVelocity<T>(parameters) );
    setBoundaryDensity (
            lattice, outlet,
            ConstantDensity<T>(1.) );
    initializeAtEquilibrium (
            lattice, lattice.getBoundingBox(),
            PoiseuilleVelocityAndDensity<T>(parameters) );

    defineDynamics(lattice, boolMask, new Porosity_GuoExternalForceBGKdynamics<T,DESCRIPTOR>(omega,rockporosity, KVC, k0), 2);
    defineDynamics(lattice, boolMask, new Porosity_GuoExternalForceBGKdynamics<T,DESCRIPTOR>(omega,cokeporosity, KVC, k0), 3);

    //pcout << "3" << endl;
    lattice.initialize();
}

void writeGif(MultiBlockLattice2D<T,DESCRIPTOR>& lattice, plint iter)
{
    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif(createFileName("u", iter, 6),
                               *computeVelocityNorm(lattice) );
}

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

int main(int argc, char* argv[]) 
{
    plbInit(&argc, &argv);
    auto outputdir = "./tmp/";
    global::directories().setOutputDir(outputdir);

    IncomprFlowParam<T> parameters(
            (T) 1e-2,  // uMax
            (T) 5,  // Re
            100,        // N
            3.39,        // lx
            1.49         // ly 
    );
    
    T KVC      = ((T)1/(T)3)* (parameters.getTau()-(T)0.5);//注意粘性系数与tau之间的关系
    T k0       = 1;//.e-12;应该单位换算之后取格子单位
    T rockporosity = 0.01;//porosity can not = 0 (invK 必须有意义)
    T cokeporosity = 0.7;

    const T logT     = (T)0.02;
    const T imSave   = (T)0.01;
    const T vtkSave  = (T)0.5;
    const T maxT     = (T)50.1;
    pcout << parameters.getTau() << endl;
    pcout << parameters.getNx() << endl;
    pcout << parameters.getNy() << endl;
    writeLogFile(parameters, "Poiseuille flow");

    MultiBlockLattice2D<T, DESCRIPTOR> lattice (parameters.getNx(), parameters.getNy(),
            new Porosity_GuoExternalForceBGKdynamics<T,DESCRIPTOR>(parameters.getOmega(),1, KVC, k0) );
    
    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();
    
    MultiScalarField2D<int> boolMask(parameters.getNx(), parameters.getNy());  
    plb_ifstream ifile("tiny_structure.txt");
    ifile >> boolMask;

    cylinderSetup(lattice, parameters, *boundaryCondition, parameters.getOmega(), rockporosity, cokeporosity, KVC, k0, boolMask);
    
    // Main loop over time iterations.
    for (plint iT=0; iT*parameters.getDeltaT()<maxT; ++iT) {
        // At this point, the state of the lattice corresponds to the
        //   discrete time iT. However, the stored averages (getStoredAverageEnergy
        //   and getStoredAverageDensity) correspond to the previous time iT-1.

       if (iT%parameters.nStep(imSave)==0) {
            pcout << "; Saving Gif ..." << endl;
            writeGif(lattice, iT);
        }

        if (iT%parameters.nStep(vtkSave)==0 && iT>0) {
            pcout << "Saving VTK file ..." << endl;
            writeVTK(lattice, parameters, iT);
        }

        if (iT%parameters.nStep(logT)==0) {
            pcout << "step " << iT
                  << "; t=" << iT*parameters.getDeltaT();
        }


        if (iT%parameters.nStep(logT)==0) {
            pcout << "step " << iT
                  << "; t=" << iT*parameters.getDeltaT()
                  << "; RMS error=" << computeRMSerror(lattice, parameters);
            }
        // Lattice Boltzmann iteration step.
        lattice.collideAndStream();

        // At this point, the state of the lattice corresponds to the
        //   discrete time iT+1, and the stored averages are upgraded to time iT.
        if (iT%parameters.nStep(logT)==0) {
            pcout << "; av energy ="
                  << setprecision(10) << getStoredAverageEnergy<T>(lattice)
                  << "; av rho ="
                  << getStoredAverageDensity<T>(lattice) << endl;
        }
    }
    
    delete boundaryCondition;
}

