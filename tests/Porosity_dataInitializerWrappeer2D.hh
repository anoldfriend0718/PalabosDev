/** \file
 * Helper functions for domain initialization -- generic implementation.
 */
#ifndef POROSITY_DATA_INITIALIZER_WRAPPER_2D_HH
#define POROSITY_DATA_INITIALIZER_WRAPPER_2D_HH


#include "dataProcessors/dataInitializerWrapper2D.h"
#include "core/cell.h"
#include "atomicBlock/dataProcessorWrapper2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataField2D.h"
#include "multiBlock/multiBlockLattice2D.h"
#include "multiBlock/multiDataField2D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"


namespace plb {

template<typename T, template<class U> class Descriptor, class Function>
void setPorosity(MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain, Function f) {
    applyIndexed(lattice, domain, new SetCustomPorosityFunctional2D<T,Descriptor,Function>(f) );
}

}  // namespace plb

#endif  // POROSITY_DATA_INITIALIZER_WRAPPER_2D_HH