/** \file
 * Helper functions for data field initialization -- generic implementation.
 */
#ifndef POROSITY_DATA_INITIALIZER_GENERICS_2D_H
#define POROSITY_DATA_INITIALIZER_GENERICS_2D_H


#include "core/globalDefs.h"
#include "atomicBlock/dataProcessorWrapper2D.h"
#include "multiBlock/multiDataProcessorWrapper2D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "multiGrid/multiGridUtil.h"


namespace plb {

    /* ************ Class SetCustomPorosityFunctional2D ********** */

template<typename T, template<typename U> class Descriptor, class PorosityFunction>
SetCustomPorosityFunctional2D<T,Descriptor,PorosityFunction>::
    SetCustomPorosityFunctional2D(PorosityFunction f_)
        : f(f_)
{ }

template<typename T, template<typename U> class Descriptor, class PorosityFunction>
void SetCustomPorosityFunctional2D<T,Descriptor,PorosityFunction>::execute (
        plint iX, plint iY, Cell<T,Descriptor>& cell ) const
{
    T porosity = f(iX, iY);
    cell.getDynamics().setPorosity(porosity);
}

template<typename T, template<typename U> class Descriptor, class PorosityFunction>
SetCustomPorosityFunctional2D<T,Descriptor,PorosityFunction>*
    SetCustomPorosityFunctional2D<T,Descriptor,PorosityFunction>::clone() const
{
    return new SetCustomPorosityFunctional2D<T,Descriptor,PorosityFunction>(*this);
}

}  // namespace plb

#endif  // POROSITY_DATA_INITIALIZER_GENERICS_2D_H
