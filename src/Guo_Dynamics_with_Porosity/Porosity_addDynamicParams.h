#ifndef POROSITY_DYNAMICS_H
#define POROSITY_DYNAMICS_H

#include "core/dynamics.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "core/blockStatistics.h"
#include "core/array.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

namespace dynamicParams {
    //porosity for Porosity_BGKdynamics
    const plint porosity = 1;
    //KVC = Kinematic viscosity coefficient
    const plint KVC = 1;
    //k0 is the coefficient of the permeability K
    const plint k0 = 1;
}

}  // namespace plb

#endif  // POROSITY_DYNAMICS_H