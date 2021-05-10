#ifndef MT_DYNAMICS_HH
#define MT_DYNAMICS_HH

#include "core/latticeStatistics.h"
#include "core/dynamicsIdentifiers.h"
#include "complexDynamics/advectionDiffusionDynamics.h"
#include "latticeBoltzmann/advectionDiffusionMomentTemplates.h"
#include "latticeBoltzmann/advectionDiffusionDynamicsTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/offEquilibriumAdvectionDiffusionTemplates.h"

#include "MTdynamics.h"
#include "MTdynamicsTemplates.h"

namespace plb {

/* *************** Class MTdynamics *********************************************** */

template<typename T, template<typename U> class Descriptor>
T MTdynamics<T,Descriptor>::getParameter(plint whichParameter) const {
    switch (whichParameter) {
        case dynamicParams::porosity   : return this->getPorosity();
    };
    return 0.;
}

template<typename T, template<typename U> class Descriptor>
void MTdynamics<T,Descriptor>::setParameter(plint whichParameter, T value) {
    switch (whichParameter) {
        case dynamicParams::porosity   : setPorosity(value);
    };
}

template<typename T, template<typename U> class Descriptor>
T MTdynamics<T,Descriptor>::getPorosity() const {
    return porosity;
}

template<typename T, template<typename U> class Descriptor>
void MTdynamics<T,Descriptor>::setPorosity(T porosity_) {
    porosity = porosity_;
}

template<typename T, template<typename U> class Descriptor>
int MTdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,MTdynamics<T,Descriptor> >("MT");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
MTdynamics<T,Descriptor>::MTdynamics(T omega_,T porosity_)
    : AdvectionDiffusionDynamics<T,Descriptor>(omega_),porosity(porosity_)
{ }

template<typename T, template<typename U> class Descriptor>
MTdynamics<T,Descriptor>::MTdynamics(HierarchicUnserializer& unserializer)
    : AdvectionDiffusionDynamics<T,Descriptor>(T()),porosity(T())
{
    this->unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
MTdynamics<T,Descriptor>* MTdynamics<T,Descriptor>::clone() const {
    return new MTdynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int MTdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void MTdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics& statistics )
{
    T rhoBar;//一般而言rhoBar=rho-1,但此处其意义为porosity*C-1
    Array<T,Descriptor<T>::d> jEq;//此处意义为rho*u,u依靠外界引入
    advectionDiffusionMomentTemplates<T,Descriptor>::get_rhoBar_jEq(cell, rhoBar, jEq);
    
    T uSqr = MTdynamicsTemplates<T,Descriptor>::
             MTcollision(cell, rhoBar, jEq, this->getOmega(), this->getPorosity());
    
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
void MTdynamics<T,Descriptor>::collideExternal (
        Cell<T,Descriptor>& cell, T rhoBar,
        Array<T,Descriptor<T>::d> const& jEq, T thetaBar, BlockStatistics& stat )
{
    T uSqr = MTdynamicsTemplates<T,Descriptor>::MTcollision(cell, rhoBar, jEq, this->getOmega(), this->getPorosity());
    
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
T MTdynamics<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& jEq,
                                                T jSqr, T thetaBar) const
{
    return MTdynamicsTemplates<T,Descriptor>::MTequilibrium(iPop, rhoBar, jEq, this->getPorosity());
}

}

#endif  // MT_DYNAMICS_HH