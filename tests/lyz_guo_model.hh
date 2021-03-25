/** \file
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- generic implementation.
 */
#ifndef lyz_guo_model_HH
#define lyz_guo_model_HH

#include "basicDynamics/lyz_guo_model.h"
#include "core/cell.h"
#include "core/dynamicsIdentifiers.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/externalForceTemplates.h"
#include "latticeBoltzmann/offEquilibriumTemplates.h"
#include "latticeBoltzmann/d3q13Templates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "core/latticeStatistics.h"
#include <algorithm>
#include <limits>
#include <cstdlib>

#include "latticeBoltzmann/NewDynamicsTemplates.h"

namespace plb {

/* *************** Class guo_BGKdynamics *********************************************** */
template<typename T, template<typename U> class Descriptor>
int guo_BGKdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,guo_BGKdynamics<T,Descriptor> >("guo_BGK");
/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
guo_BGKdynamics<T,Descriptor>::guo_BGKdynamics(T omega_ )
    : IsoThermalBulkDynamics<T,Descriptor>(omega_)
{ }
template<typename T, template<typename U> class Descriptor>
guo_BGKdynamics<T,Descriptor>::guo_BGKdynamics(HierarchicUnserializer& unserializer)
    : IsoThermalBulkDynamics<T,Descriptor>(T())
{
    this->unserialize(unserializer);
}
template<typename T, template<typename U> class Descriptor>
guo_BGKdynamics<T,Descriptor>* guo_BGKdynamics<T,Descriptor>::clone() const {
    return new guo_BGKdynamics<T,Descriptor>(*this);
}
template<typename T, template<typename U> class Descriptor>
int guo_BGKdynamics<T,Descriptor>::getId() const {
    return id;
}
template<typename T, template<typename U> class Descriptor>
void guo_BGKdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics& statistics )
{
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T uSqr = dynamicsTemplates<T,Descriptor>::guo_bgk_collision(cell, rhoBar, j, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}
template<typename T, template<typename U> class Descriptor>
void guo_BGKdynamics<T,Descriptor>::collideExternal (
        Cell<T,Descriptor>& cell, T rhoBar,
        Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat )
{
    T uSqr = dynamicsTemplates<T,Descriptor>::guo_bgk_collision(cell, rhoBar, j, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}
template<typename T, template<typename U> class Descriptor>
T guo_BGKdynamics<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                                T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::guo_bgk_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}
template<typename T, template<typename U> class Descriptor>
void guo_BGKdynamics<T,Descriptor>::decomposeOrder0 (
        Cell<T,Descriptor> const& cell, std::vector<T>& rawData ) const
{
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);
    rawData[0] = rhoBar;
    j.to_cArray(&rawData[1]);
    
    Array<T,Descriptor<T>::q> fEq;
    dynamicsTemplates<T,Descriptor>::guo_bgk_equilibria( rhoBar, Descriptor<T>::invRho(rhoBar), j, jSqr, fEq );
    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
        rawData[1+Descriptor<T>::d+iPop] =
            cell[iPop] - fEq[iPop];
    }
    int offset = 1+Descriptor<T>::d+Descriptor<T>::q;
    for (plint iExt=0; iExt<Descriptor<T>::ExternalField::numScalars; ++iExt) {
        rawData[offset+iExt] = *cell.getExternal(iExt);
    }
}
template<typename T, template<typename U> class Descriptor>
void guo_BGKdynamics<T,Descriptor>::recomposeOrder0 (
        Cell<T,Descriptor>& cell, std::vector<T> const& rawData ) const
{
    T rhoBar = rawData[0];
    Array<T,Descriptor<T>::d> j;
    j.from_cArray(&rawData[1]);
    T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);
    
    Array<T,Descriptor<T>::q> fEq;
    dynamicsTemplates<T,Descriptor>::guo_bgk_equilibria( rhoBar, Descriptor<T>::invRho(rhoBar), j, jSqr, fEq );
    
    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
        cell[iPop] = fEq[iPop] + rawData[1+Descriptor<T>::d+iPop];
    }
    int offset = 1+Descriptor<T>::d+Descriptor<T>::q;
    for (plint iExt=0; iExt<Descriptor<T>::ExternalField::numScalars; ++iExt) {
        *cell.getExternal(iExt) = rawData[offset+iExt];
    }
}


}  // namespace plb

#endif  // ISO_THERMAL_DYNAMICS_HH

