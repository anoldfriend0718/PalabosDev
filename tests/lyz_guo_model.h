/** \file
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- header file.
 */
#ifndef lyz_guo_model_H
#define lyz_guo_model_H

#include "core/globalDefs.h"
#include "core/dynamics.h"

namespace plb {

/// Implementation of guo_BGK dynamics
template<typename T, template<typename U> class Descriptor>
class guo_BGKdynamics : public IsoThermalBulkDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    guo_BGKdynamics(T omega_,T eps);
    guo_BGKdynamics(HierarchicUnserializer& unserializer);

    /// Clone the object on its dynamic type.
    virtual guo_BGKdynamics<T,Descriptor>* clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    T getPorosity() {return eps;}

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
private:
    virtual void decomposeOrder0(Cell<T,Descriptor> const& cell, std::vector<T>& rawData) const;
    virtual void recomposeOrder0(Cell<T,Descriptor>& cell, std::vector<T> const& rawData) const;
private:
    static int id;

    T eps; //porosity
};



}  // namespace plb

#endif  // ISO_THERMAL_DYNAMICS_H
