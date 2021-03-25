

 #include "palabos2D.h"
 #include "palabos2D.hh"
 //#include "basicDynamics/lyz_guo_model.h"
 //#include "basicDynamics/lyz_guo_model.hh"
 //#include "latticeBoltzmann/NewDynamicsTemplates.h"
 #include <iostream>
 #include <iomanip>

 using namespace plb;
 using namespace std;


//13 
typedef double T;
//14
 #define DESCRIPTOR plb::descriptors::D2Q9Descriptor

void defineInitialDensityAtCenter(MultiBlockLattice2D<T,DESCRIPTOR>& lattice)
 {
 // The lattice is of size nx-by-ny
 const plint nx = lattice.getNx();
 const plint ny = lattice.getNy();

 // Create a Box2D which describes the location of cells with a slightly
 // higher density.
 plint centralSquareRadius = nx/6;
 plint centerX = nx/3;
 plint centerY = ny/4;
 Box2D centralSquare (
 centerX - centralSquareRadius, centerX + centralSquareRadius,
 centerY - centralSquareRadius, centerY + centralSquareRadius );

 // All cells have initially density rho ...
 T rho0 = 1.;
 // .. except for those in the box "centralSquare" which have density
 // rho+deltaRho
 T deltaRho = 1.e-4;
 Array<T,2> u0(0,0);

 // Initialize constant density everywhere.
 initializeAtEquilibrium (
 lattice, lattice.getBoundingBox(), rho0, u0 );
 // And slightly higher density in the central box.
 initializeAtEquilibrium (
 lattice, centralSquare, rho0 + deltaRho, u0 );

 lattice.initialize();
 }

 int main(int argc, char* argv[]) {
 plbInit(&argc, &argv);
 global::directories().setOutputDir("./tmp_lyz_guo/");

 const plint maxIter = 1000; // Iterate during 1000 steps.
 const plint nx = 600; // Choice of lattice dimensions.
 const plint ny = 600;
 const T omega = 1.; // Choice of the relaxation parameter

 MultiBlockLattice2D<T, DESCRIPTOR> lattice (
 nx, ny, new CompleteRegularizedBGKdynamics<T,DESCRIPTOR>(omega) );

 lattice.periodicity().toggleAll(true); // Use periodic boundaries.

 defineInitialDensityAtCenter(lattice);

 // Main loop over time iterations.
 for (plint iT=0; iT<maxIter; ++iT) {
 if (iT%40==0) { // Write an image every 40th time step.
 pcout << "Writing GIF file at iT=" << iT << endl;
 // Instantiate an image writer with the color map "leeloo".
 ImageWriter<T> imageWriter("leeloo");
 // Write a GIF file with colors rescaled to the range of values
 // in the matrix
 imageWriter.writeScaledGif (
 createFileName("u", iT, 6),
 *computeVelocityNorm(lattice) );
 }
 // Execute lattice Boltzmann iteration.
 lattice.collideAndStream();
 }
 }