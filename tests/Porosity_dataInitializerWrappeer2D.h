/** \file
 * Helper functions for domain initialization -- header file.
 */
#ifndef POROSITY_DATA_INITIALIZER_WRAPPER_2D_H
#define POROSITY_DATA_INITIALIZER_WRAPPER_2D_H


#include "core/globalDefs.h"
#include "atomicBlock/dataProcessorWrapper2D.h"
#include "multiBlock/multiDataProcessorWrapper2D.h"
#include "core/dynamics.h"
#include "dataProcessors/dataInitializerFunctional2D.h"
#include "sitmo/prng_engine.hpp"


namespace plb {

template<typename T, template<class U> class Descriptor, class Function>
void setPorosity(MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain, Function f);


}  // namespace plb

#endif  // POROSITY_DATA_INITIALIZER_WRAPPER_2D_H