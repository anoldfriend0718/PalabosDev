/** \file
 * 2D specialization of dynamicsTemplates functions.
 */

#ifndef NEW_DYNAMICS_TEMPLATES_2D_H
#define NEW_DYNAMICS_TEMPLATES_2D_H

#include "core/globalDefs.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {



// Efficient specialization for D2Q9 base lattice
template<typename T> struct NewdynamicsTemplatesImpl<T, descriptors::D2Q9DescriptorBase<T> > {

typedef descriptors::D2Q9DescriptorBase<T> D;



static T guo_bgk_equilibrium(plint iPop, T rhoBar, T invRho, Array<T,2> const& j, T jSqr, T eip ) {
    typedef descriptors::D2Q9DescriptorBase<T> L;
    T c_j = L::c[iPop][0]*j[0] + L::c[iPop][1]*j[1];
    return L::t[iPop] * ( rhoBar +
               (T)3*c_j + invRho/eip*(4.5*c_j*c_j - 1.5*jSqr) );
}



static void guo_bgk_equilibria( T rhoBar, T invRho, Array<T,D::d> const& j,
                                T jSqr, Array<T,D::q>& eqPop, T eip )
{
    T t0 = D::t[0];
    T t1 = D::t[1];
    T t2 = D::t[2];

    T kx     = (T)3 * j[0];
    T ky     = (T)3 * j[1];
    T kxSqr_ = invRho / ((T)2*eip) * kx*kx;
    T kySqr_ = invRho / ((T)2*eip) * ky*ky;
    T kxky_  = invRho /eip * kx*ky;

    T C1 = rhoBar + invRho/eip *(T)3*jSqr;
    T C2, C3;

    // i=0
    C3 = -kxSqr_ - kySqr_;
    eqPop[0] = t0 * (C1+C3);

    // i=1 and i=5
    C2 = -kx + ky;
    C3 = -kxky_;
    eqPop[1] = t1 * (C1+C2+C3);
    eqPop[5] = t1 * (C1-C2+C3);

    // i=2 and i=6
    C2 = -kx;
    C3 = -kySqr_;
    eqPop[2] = t2 * (C1+C2+C3);
    eqPop[6] = t2 * (C1-C2+C3);

    // i=3 and i=7
    C2 = -kx - ky;
    C3 = kxky_;
    eqPop[3] = t1 * (C1+C2+C3);
    eqPop[7] = t1 * (C1-C2+C3);

    // i=4 and i=8
    C2 = -ky;
    C3 = -kxSqr_;
    eqPop[4] = t2 * (C1+C2+C3);
    eqPop[8] = t2 * (C1-C2+C3);
}






static T guo_bgk_collision_base(Array<T,D::q>& f, T rhoBar, Array<T,2> const& j, T omega, T invRho, T eip ) {
    T one_m_omega = (T)1 - omega;
    T t0_omega = D::t[0] * omega;
    T t1_omega = D::t[1] * omega;
    T t2_omega = D::t[2] * omega;

    T jSqr   = j[0]*j[0] + j[1]*j[1];
    T kx     = (T)3 * j[0];
    T ky     = (T)3 * j[1];
    T kxSqr_ = invRho / ((T)2*eip) * kx*kx;
    T kySqr_ = invRho / ((T)2*eip) * ky*ky;
    T kxky_  = invRho / eip * kx*ky;

    T C1 = rhoBar + invRho/eip*(T)3*jSqr;
    T C2, C3;

    // i=0
    C3 = -kxSqr_ - kySqr_;
    f[0] *= one_m_omega; f[0] += t0_omega * (C1+C3);

    // i=1 and i=5
    C2 = -kx + ky;
    C3 = -kxky_;
    f[1] *= one_m_omega; f[1] += t1_omega * (C1+C2+C3);
    f[5] *= one_m_omega; f[5] += t1_omega * (C1-C2+C3);

    // i=2 and i=6
    C2 = -kx;
    C3 = -kySqr_;
    f[2] *= one_m_omega; f[2] += t2_omega * (C1+C2+C3);
    f[6] *= one_m_omega; f[6] += t2_omega * (C1-C2+C3);

    // i=3 and i=7
    C2 = -kx - ky;
    C3 = kxky_;
    f[3] *= one_m_omega; f[3] += t1_omega * (C1+C2+C3);
    f[7] *= one_m_omega; f[7] += t1_omega * (C1-C2+C3);

    // i=4 and i=8
    C2 = -ky;
    C3 = -kxSqr_;
    f[4] *= one_m_omega; f[4] += t2_omega * (C1+C2+C3);
    f[8] *= one_m_omega; f[8] += t2_omega * (C1-C2+C3);

    return invRho*invRho*jSqr;
}




static T guo_bgk_collision(Array<T,D::q>& f, T rhoBar, Array<T,2> const& j, T omega, T eip) {
    return guo_bgk_collision_base(f, rhoBar, j, omega, D::invRho(rhoBar), eip);
}






};  //struct dynamicsTemplatesImpl<D2Q9DescriptorBase>

}  // namespace plb

#endif  // NEW_DYNAMICS_TEMPLATES_2D_H
