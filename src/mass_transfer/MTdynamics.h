#ifndef MT_DYNAMICS_H
#define MT_DYNAMICS_H

#include "core/globalDefs.h"
#include "core/dynamics.h"

#include "Guo_Dynamics_with_Porosity/Porosity_addDynamicParams.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
class MTdynamics : public AdvectionDiffusionDynamics <T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    MTdynamics(T omega_,T porosity_);
    MTdynamics(HierarchicUnserializer& unserializer);

    /// Clone the object on its dynamic type.
    virtual MTdynamics<T,Descriptor>* clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

/* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics& statistics_);
     /// Implementation of the collision step, with imposed macroscopic variables
    /// The arguments:
    /// - rhoBar: the "rhoBar" version of the scalar rho.
    /// - jEq: the equilibrium part of the second-order moment. jEq = u*rho, where u is the external convective term.

    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(Cell<T,Descriptor>& cell, T rhoBar,
                         Array<T,Descriptor<T>::d> const& jEq, T thetaBar, BlockStatistics& stat);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& jEq,
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
}  // namespace plb

#endif  // MT_DYNAMICS_H