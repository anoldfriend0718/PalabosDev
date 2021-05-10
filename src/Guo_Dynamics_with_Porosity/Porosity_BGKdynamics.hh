#ifndef POROSITY_BGKDYNAMICS_HH
#define POROSITY_BGKDYNAMICS_HH

#include "core/globalDefs.h"
#include "core/dynamics.h"
#include "core/util.h"
#include "core/blockStatistics.h"
#include "core/array.h"
#include "basicDynamics/isoThermalDynamics.h"
#include "basicDynamics/isoThermalDynamics.hh"
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
#include "core/hierarchicSerializer.h"
#include "latticeBoltzmann/indexTemplates.h"
#include "multiGrid/multiGridUtil.h"

#include "Porosity_BGKdynamics.h"
#include "Porosity_dynamicsTemplates.h" //XQH comment: test is a bad name
#include "Porosity_dynamicsTemplates2D.h"
#include "Porosity_addDynamicParams.h"

namespace plb{

/* *************** Class Porosity_BGKdynamics ********************************** */
template<typename T, template<typename U> class Descriptor>
int Porosity_BGKdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,Porosity_BGKdynamics<T,Descriptor> >("Porosity_BGK");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
Porosity_BGKdynamics<T,Descriptor>::Porosity_BGKdynamics(T omega_,T porosity_)
    : IsoThermalBulkDynamics<T,Descriptor>(omega_),porosity(porosity_)
{ }

template<typename T, template<typename U> class Descriptor>
Porosity_BGKdynamics<T,Descriptor>::Porosity_BGKdynamics(HierarchicUnserializer& unserializer)
    : IsoThermalBulkDynamics<T,Descriptor>(T()),porosity(T())
{
    this->unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
void Porosity_BGKdynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    IsoThermalBulkDynamics<T,Descriptor>::serialize(serializer);
    serializer.addValue(porosity);
}

template<typename T, template<typename U> class Descriptor>
void Porosity_BGKdynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    IsoThermalBulkDynamics<T,Descriptor>::unserialize(unserializer);
    porosity = unserializer.readValue<T>();
    unserializer.readValue(porosity);
}

template<typename T, template<typename U> class Descriptor>
Porosity_BGKdynamics<T,Descriptor>* Porosity_BGKdynamics<T,Descriptor>::clone() const {
    return new Porosity_BGKdynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int Porosity_BGKdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void Porosity_BGKdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics& statistics )
{
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, j);

    T uSqr = Porosity_dynamicsTemplates<T,Descriptor>::
    Porosity_bgk_collision(cell, rhoBar, j, this->getOmega(), this->getPorosity());

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
void Porosity_BGKdynamics<T,Descriptor>::collideExternal (
        Cell<T,Descriptor>& cell, T rhoBar,
        Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat )
{
    T uSqr = Porosity_dynamicsTemplates<T,Descriptor>::Porosity_bgk_collision(cell, rhoBar, j, this->getOmega(), this->getPorosity());
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
T Porosity_BGKdynamics<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                                T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return Porosity_dynamicsTemplates<T,Descriptor>::Porosity_bgk_equilibrium(iPop, rhoBar, invRho, j, jSqr, this->getPorosity());
}


template<typename T, template<typename U> class Descriptor>
T Porosity_BGKdynamics<T,Descriptor>::getParameter(plint whichParameter) const {
    switch (whichParameter) {
        case dynamicParams::porosity   : return this->getPorosity();
    };
    return 0.;
}

template<typename T, template<typename U> class Descriptor>
void Porosity_BGKdynamics<T,Descriptor>::setParameter(plint whichParameter, T value) {
    switch (whichParameter) {
        case dynamicParams::porosity   : setPorosity(value);
    };
}

template<typename T, template<typename U> class Descriptor>
T Porosity_BGKdynamics<T,Descriptor>::getPorosity() const {
    return porosity;
}

template<typename T, template<typename U> class Descriptor>
void Porosity_BGKdynamics<T,Descriptor>::setPorosity(T porosity_) {
    porosity = porosity_;
}


}  // namespace plb

#endif  // POROSITY_BGKDYNAMICS_H