/** \file
 * Functionals for domain initialization -- header file.
 */
#ifndef POROSITY_DATA_INITIALIZER_FUNCTIONAL_2D_H
#define POROSITY_DATA_INITIALIZER_FUNCTIONAL_2D_H


#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "core/dynamics.h"
#include "sitmo/prng_engine.hpp"




namespace plb {

/* ************* Class SetCustomPorosityFunctional2D ******************* */

template<typename T, template<typename U> class Descriptor, class PorosityFunction>
class SetCustomPorosityFunctional2D : public OneCellIndexedFunctional2D<T,Descriptor>
{
public:
    SetCustomPorosityFunctional2D(PorosityFunction f_);
    virtual void execute(plint iX, plint iY, Cell<T,Descriptor>& cell) const;
    virtual SetCustomPorosityFunctional2D<T,Descriptor,PorosityFunction>* clone() const;
private:
    PorosityFunction f;
};

}  // namespace plb

#endif  // POROSITY_DATA_INITIALIZER_FUNCTIONAL_2D_H

// Explicitly include generic algorithms which are never precompiled (not even in precompiled version)
/*
#include "dataProcessors/dataInitializerGenerics2D.h"
*/
