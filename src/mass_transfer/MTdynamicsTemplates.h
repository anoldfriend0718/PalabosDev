#ifndef MT_DYNAMICS_TEMPLATES_H
#define MT_DYNAMICS_TEMPLATES_H

#include "core/globalDefs.h"
#include "core/cell.h"
#include "core/util.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/offEquilibriumTemplates.h"
#include "latticeBoltzmann/offEquilibriumAdvectionDiffusionTemplates.h"

#include "Guo_Dynamics_with_Porosity/Porosity_addDynamicParams.h"

namespace plb {

template<typename T, class Descriptor> struct MTdynamicsTemplatesImpl;

/// This structure forwards the calls to the appropriate helper class
template<typename T, template<typename U> class Descriptor>
struct MTdynamicsTemplates {

static T MTequilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& jEq, T porosity) 
{
    return MTdynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::MTequilibrium(iPop, rhoBar, jEq, porosity);
}

static T MTcollision(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& jEq, T omega, T porosity)
{
    return MTdynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::MTcollision(cell.getRawPopulations(), rhoBar, jEq, omega, porosity);
}

};  // struct MTdynamicsTemplates


/// All helper functions are inside this structure
template<typename T, class Descriptor>
struct MTdynamicsTemplatesImpl {

static T MTequilibrium(plint iPop, T rhoBar, Array<T,Descriptor::d> const& jEq, T porosity) 
{
    T c_j = Descriptor::c[iPop][0]*jEq[0];
    for (int iD=1; iD < Descriptor::d; ++iD) {
       c_j += Descriptor::c[iPop][iD]*jEq[iD];
    }
    return Descriptor::t[iPop] * (rhoBar + Descriptor::invCs2 * c_j / porosity);
}

static T MTcollision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& jEq, T omega, T porosity) 
{
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = VectorTemplateImpl<T,Descriptor::d>::normSqr(jEq);
    for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
        f[iPop] *= (T)1-omega;
        f[iPop] += omega * MTdynamicsTemplatesImpl<T,Descriptor>::MTequilibrium (
                                iPop, rhoBar, jEq, porosity );
    }
    return jSqr*invRho*invRho;
}


};  // struct MTdynamicsTemplatesImpl

}  // namespace plb

#include "MTdynamicsTemplates2D.h"


#endif  // MT_DYNAMICS_TEMPLATES_H