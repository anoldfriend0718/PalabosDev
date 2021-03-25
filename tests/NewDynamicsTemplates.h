/** \file
 * Helper functions for the implementation of LB dynamics. This file is all
 * about efficiency. The generic template code is specialized for commonly
 * used Lattices, so that a maximum performance can be taken out of each
 * case.
 */
#ifndef NEW_DYNAMICS_TEMPLATES_H
#define NEW_DYNAMICS_TEMPLATES_H

#include "core/globalDefs.h"
#include "core/cell.h"
#include "core/util.h"
#include "latticeBoltzmann/offEquilibriumTemplates.h"

namespace plb {

template<typename T, class Descriptor> struct dynamicsTemplatesImpl;

/// This structure forwards the calls to the appropriate helper class
template<typename T, template<typename U> class Descriptor>
struct dynamicsTemplates {

static T guo_bgk_equilibrium(plint iPop, T rhoBar, T invRho, Array<T,Descriptor<T>::d> const& j, T jSqr) {
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::guo_bgk_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}
static void guo_bgk_equilibria( T rhoBar, T invRho, Array<T,Descriptor<T>::d> const& j,
                                T jSqr, Array<T,Descriptor<T>::q>& eqPop )
{
    dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::guo_bgk_equilibria(rhoBar, invRho, j, jSqr, eqPop);
}
static T guo_bgk_collision(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j, T omega)
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::guo_bgk_collision(cell.getRawPopulations(), rhoBar, j, omega);
}





};  // struct dynamicsTemplates


/// All helper functions are inside this structure
template<typename T, class Descriptor>
struct dynamicsTemplatesImpl {
//eip for 孔隙率
static T guo_bgk_equilibrium(plint iPop, T rhoBar, T invRho, Array<T,Descriptor::d> const& j, T jSqr) {
    T c_j = Descriptor::c[iPop][0]*j[0];
    for (int iD=1; iD < Descriptor::d; ++iD) {
       c_j += Descriptor::c[iPop][iD]*j[iD];
    }
    return Descriptor::t[iPop] *(rhoBar+(T)1) (
           (T)1 + Descriptor::invCs2 * c_j +
           Descriptor::invCs2/((T)2*eip) * (
               Descriptor::invCs2 * c_j*c_j - jSqr )
       );
}
static void guo_bgk_ma2_equilibria( T rhoBar, T invRho, Array<T,Descriptor::d>
        const& j,
                                T jSqr, Array<T,Descriptor::q>& eqPop )
{
    for (int iPop=0; iPop<Descriptor::q; ++iPop) {
        eqPop[iPop] = guo_bgk_equilibrium(iPop, rhoBar, invRho, j, jSqr);
    }
}
static T guo_bgk_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& j, T omega) {
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = VectorTemplateImpl<T,Descriptor::d>::normSqr(j);
    for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
        f[iPop] *= (T)1-omega;
        f[iPop] += omega * dynamicsTemplatesImpl<T,Descriptor>::guo_bgk_equilibrium (
                                iPop, rhoBar, invRho, j, jSqr );
    }
    return jSqr*invRho*invRho;
}


};  // struct dynamicsTemplatesImpl

}  // namespace plb

#include "latticeBoltzmann/dynamicsTemplates2D.h"
#include "latticeBoltzmann/dynamicsTemplates3D.h"

#endif  // DYNAMICS_TEMPLATES_H

