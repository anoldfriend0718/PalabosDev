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
}

}  // namespace plb

#endif  // POROSITY_DYNAMICS_H