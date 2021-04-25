#ifndef POROSITY_DYNAMICS_TEMPLATES_H
#define POROSITY_DYNAMICS_TEMPLATES_H

#include "core/globalDefs.h"
#include "core/cell.h"
#include "core/util.h"
#include "latticeBoltzmann/offEquilibriumTemplates.h"
#include "test/Porosity_dynamics.h"

namespace plb {

template<typename T, class Descriptor> struct Porosity_dynamicsTemplatesImpl;

/// This structure forwards the calls to the appropriate helper class
template<typename T, template<typename U> class Descriptor>
struct Porosity_dynamicsTemplates {

static T Porosity_bgk_equilibrium(plint iPop, T rhoBar, T invRho, Array<T,Descriptor<T>::d> const& j, T jSqr, T porosity) {
    return Porosity_dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::Porosity_bgk_equilibrium(iPop, rhoBar, invRho, j, jSqr, porosity);
}

static void Porosity_bgk_equilibria( T rhoBar, T invRho, Array<T,Descriptor<T>::d> const& j,
                                T jSqr, Array<T,Descriptor<T>::q>& eqPop, T porosity )
{
    Porosity_dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::Porosity_bgk_equilibria(rhoBar, invRho, j, jSqr, eqPop, porosity);
}

static T Porosity_bgk_collision(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j, T omega, T porosity)
{
    return Porosity_dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::Porosity_bgk_collision(cell.getRawPopulations(), rhoBar, j, omega, porosity);
}

};  // struct Porosity_dynamicsTemplates


/// All helper functions are inside this structure
template<typename T, class Descriptor>
struct Porosity_dynamicsTemplatesImpl {

static T Porosity_bgk_equilibrium(plint iPop, T rhoBar, T invRho, Array<T,Descriptor::d> const& j, T jSqr, T porosity) {
    T c_j = Descriptor::c[iPop][0]*j[0];
    for (int iD=1; iD < Descriptor::d; ++iD) {
       c_j += Descriptor::c[iPop][iD]*j[iD];
    }
    return Descriptor::t[iPop] * (
           rhoBar + Descriptor::invCs2 * c_j +
           Descriptor::invCs2/((T)2 * porosity) * invRho * (
               Descriptor::invCs2 * c_j*c_j - jSqr )
       );
}

static void Porosity_bgk_equilibria( T rhoBar, T invRho, Array<T,Descriptor::d>
        const& j, T jSqr, Array<T,Descriptor::q>& eqPop, T porosity )
{
    for (int iPop=0; iPop<Descriptor::q; ++iPop) {
        eqPop[iPop] = Porosity_bgk_equilibrium(iPop, rhoBar, invRho, j, jSqr, porosity);
    }
}

static T Porosity_bgk_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& j, T omega, T porosity) {
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = VectorTemplateImpl<T,Descriptor::d>::normSqr(j);
    for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
        f[iPop] *= (T)1-omega;
        f[iPop] += omega * dynamicsTemplatesImpl<T,Descriptor>::Porosity_bgk_equilibrium (  //XQH comment :if the space name is wronge?
                                iPop, rhoBar, invRho, j, jSqr, porosity );
    }
    return jSqr*invRho*invRho;
}


};  // struct Porosity_dynamicsTemplatesImpl

}  // namespace plb

#include "test/Porosity_dynamicsTemplates2D.h"

#endif  // POROSITY_DYNAMICS_TEMPLATES_H