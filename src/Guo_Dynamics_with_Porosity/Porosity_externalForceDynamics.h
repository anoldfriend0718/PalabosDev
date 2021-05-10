#ifndef POROSITY_EXTERNAL_FORCE_DYNAMICS_H
#define POROSITY_EXTERNAL_FORCE_DYNAMICS_H

#include "core/globalDefs.h"
#include "core/dynamics.h"
#include "basicDynamics/isoThermalDynamics.h"

#include "basicDynamics/externalForceDynamics.h"
#include "basicDynamics/externalForceDynamics.hh"

#include "Porosity_addDynamicParams.h"

namespace plb {


/// Implementation of O(Ma^2) BGK dynamics with an external force (Porosity_Guo approach)
/// KVC = Kinematic viscosity coefficient
template<typename T, template<typename U> class Descriptor>
class Porosity_GuoExternalForceBGKdynamics : public ExternalForceDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    Porosity_GuoExternalForceBGKdynamics(T omega_, T porosity_, T KVC_, T k0_);
    Porosity_GuoExternalForceBGKdynamics(HierarchicUnserializer& unserializer);

    /// Clone the object on its dynamic type.
    virtual Porosity_GuoExternalForceBGKdynamics<T,Descriptor>* clone() const;

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

    T getKVC() const;
    
    void setKVC(T KVC_);

    T getk0() const;

    void setk0(T k0_);

    virtual void computeVelocity (Cell<T,Descriptor> const& cell, Array<T,Descriptor<T>::d>& u) const;

    virtual void computeVelocityExternal (Cell<T,Descriptor> const& cell, T rhoBar,
                     Array<T,Descriptor<T>::d> const& j,Array<T,Descriptor<T>::d>& u ) const;

private:
    static int id;
private:
    T porosity;
    T KVC;
    T k0;
};

}  // namespace plb

#endif  // POROSITY_EXTERNAL_FORCE_DYNAMICS_H
