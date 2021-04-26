#ifndef POROSITY_BGKDYNAMICS_H
#define POROSITY_BGKDYNAMICS_H

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

#include "test/Porosity_dynamicsTemplates.h" //XQH comment: test is a bad name
#include "test/Porosity_dynamicsTemplates2D.h"
#include "test/Porosity_dynamics.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
class Porosity_BGKdynamics : public IsoThermalBulkDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    Porosity_BGKdynamics(T omega_,T porosity_);
    Porosity_BGKdynamics(HierarchicUnserializer& unserializer);

    /// Clone the object on its dynamic type.
    virtual Porosity_BGKdynamics<T,Descriptor>* clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer& serializer) const;

    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer& unserializer);

/* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics& statistics_);

    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(Cell<T,Descriptor>& cell, T rhoBar,
                         Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
    
    /// Set local value of any generic parameter
    virtual void setParameter(plint whichParameter, T value);

    /// Get local value of any generic parameter
    virtual T getParameter(plint whichParameter) const;   

    T getPorosity() const;

    void setPorosity(T porosity_);
    
private:
    virtual void decomposeOrder0(Cell<T,Descriptor> const& cell, std::vector<T>& rawData) const;
    virtual void recomposeOrder0(Cell<T,Descriptor>& cell, std::vector<T> const& rawData) const;
private:
    static int id;
private:
    T porosity;
};

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

/* *************** Class Porosity_BGKdynamics *********************************************** */

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
void Porosity_BGKdynamics<T,Descriptor>::decomposeOrder0 (
        Cell<T,Descriptor> const& cell, std::vector<T>& rawData ) const
{
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);
    rawData[0] = rhoBar;
    j.to_cArray(&rawData[1]);
    
    Array<T,Descriptor<T>::q> fEq;
    Porosity_dynamicsTemplates<T,Descriptor>::Porosity_bgk_equilibria( rhoBar, Descriptor<T>::invRho(rhoBar), j, jSqr, fEq, this->getPorosity() );

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
void Porosity_BGKdynamics<T,Descriptor>::recomposeOrder0 (
        Cell<T,Descriptor>& cell, std::vector<T> const& rawData ) const
{
    T rhoBar = rawData[0];
    Array<T,Descriptor<T>::d> j;
    j.from_cArray(&rawData[1]);
    T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);

    
    Array<T,Descriptor<T>::q> fEq;
    Porosity_dynamicsTemplates<T,Descriptor>::Porosity_bgk_equilibria( rhoBar, Descriptor<T>::invRho(rhoBar), j, jSqr, fEq, this->getPorosity() );
    
    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
        cell[iPop] = fEq[iPop] + rawData[1+Descriptor<T>::d+iPop];
    }

    int offset = 1+Descriptor<T>::d+Descriptor<T>::q;
    for (plint iExt=0; iExt<Descriptor<T>::ExternalField::numScalars; ++iExt) {
        *cell.getExternal(iExt) = rawData[offset+iExt];
    }
}


}  // namespace plb

#endif  // POROSITY_BGKDYNAMICS_H