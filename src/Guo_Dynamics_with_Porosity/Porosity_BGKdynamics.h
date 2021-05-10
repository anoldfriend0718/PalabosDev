#ifndef POROSITY_BGKDYNAMICS_H
#define POROSITY_BGKDYNAMICS_H

#include "core/globalDefs.h"
#include "core/dynamics.h"

#include "basicDynamics/isoThermalDynamics.h"

#include "Porosity_addDynamicParams.h"

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
    static int id;
private:
    T porosity;

};

} //namespace plb

#endif 
