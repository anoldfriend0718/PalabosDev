static T bgk_ma2_collision_base(Array<T,D::q>& f, T rhoBar, Array<T,2> const& j, T omega, T invRho ) {
    T one_m_omega = (T)1 - omega;
    T t0_omega = D::t[0] * omega;
    T t1_omega = D::t[1] * omega;
    T t2_omega = D::t[2] * omega;

    T jSqr   = j[0]*j[0] + j[1]*j[1];
    T kx     = (T)3 * j[0];
    T ky     = (T)3 * j[1];
    T kxSqr_ = invRho / (T)2 * kx*kx;
    T kySqr_ = invRho / (T)2 * ky*ky;
    T kxky_  = invRho * kx*ky;

    T C1 = rhoBar + invRho*(T)3*jSqr;
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


static T complete_bgk_ma2_collision_base(Array<T,D::q>& f, T rhoBar, T invRho, Array<T,2> const& j, T omega) {
    T one_m_omega = (T)1 - omega;
    T t0_omega = D::t[0] * omega;
    T t1_omega = D::t[1] * omega;
    T t2_omega = D::t[2] * omega;

    T jSqr   = j[0]*j[0] + j[1]*j[1];
    T kx     = (T)3 * j[0];
    T ky     = (T)3 * j[1];
    T kxSqr_ = invRho / (T)2 * kx*kx;
    T kySqr_ = invRho / (T)2 * ky*ky;
    T kxky_  = invRho * kx*ky;

    T C1 = rhoBar + invRho*(T)3*jSqr;
    T C2, C3;
    
    T ux = j[0]*invRho;
    T uy = j[1]*invRho;
    
    T ux2 = ux*ux; T uy2 = uy*uy;

    // i=0
    C3 = -kxSqr_ - kySqr_;
    f[0] *= one_m_omega; f[0] += t0_omega * (C1+C3)+omega*j[0]*ux*uy2;
//     f[0] *= one_m_omega; f[0] += t0_omega * (C1+C3+j[0]*ux*uy2);

    // i=1 and i=5
    C2 = -kx + ky;
    C3 = -kxky_;
    f[1] *= one_m_omega; f[1] += t1_omega * (C1+C2+C3) + omega*(T)0.25*j[0]*uy*(ux*uy+ux-uy);
    f[5] *= one_m_omega; f[5] += t1_omega * (C1-C2+C3) + omega*(T)0.25*j[0]*uy*(ux*uy-ux+uy);
//     f[1] *= one_m_omega; f[1] += t1_omega * (C1+C2+C3+(T)0.25*j[0]*uy*(ux*uy+ux-uy));
//     f[5] *= one_m_omega; f[5] += t1_omega * (C1-C2+C3+(T)0.25*j[0]*uy*(ux*uy-ux+uy));

    // i=2 and i=6
    C2 = -kx;
    C3 = -kySqr_;
    f[2] *= one_m_omega; f[2] += t2_omega * (C1+C2+C3) - omega*(T)0.5*j[0]*uy2*(ux-(T)1);
    f[6] *= one_m_omega; f[6] += t2_omega * (C1-C2+C3) - omega*(T)0.5*j[0]*uy2*(ux+(T)1);
//     f[2] *= one_m_omega; f[2] += t2_omega * (C1+C2+C3 - (T)0.5*j[0]*uy2*(ux-(T)1));
//     f[6] *= one_m_omega; f[6] += t2_omega * (C1-C2+C3 - (T)0.5*j[0]*uy2*(ux+(T)1));

    // i=3 and i=7
    C2 = -kx - ky;
    C3 = kxky_;
    f[3] *= one_m_omega; f[3] += t1_omega * (C1+C2+C3) + omega*(T)0.25*j[0]*uy*(ux*uy-ux-uy);
    f[7] *= one_m_omega; f[7] += t1_omega * (C1-C2+C3) + omega*(T)0.25*j[0]*uy*(ux*uy+ux+uy);
//     f[3] *= one_m_omega; f[3] += t1_omega * (C1+C2+C3 + (T)0.25*j[0]*uy*(ux*uy-ux-uy));
//     f[7] *= one_m_omega; f[7] += t1_omega * (C1-C2+C3 + (T)0.25*j[0]*uy*(ux*uy+ux+uy));

    // i=4 and i=8
    C2 = -ky;
    C3 = -kxSqr_;
    f[4] *= one_m_omega; f[4] += t2_omega * (C1+C2+C3) - omega*(T)0.5*j[1]*ux2*(uy-(T)1);
    f[8] *= one_m_omega; f[8] += t2_omega * (C1-C2+C3) - omega*(T)0.5*j[1]*ux2*(uy+(T)1);
//     f[4] *= one_m_omega; f[4] += t2_omega * (C1+C2+C3 - (T)0.5*j[1]*ux2*(uy-(T)1));
//     f[8] *= one_m_omega; f[8] += t2_omega * (C1-C2+C3 - (T)0.5*j[1]*ux2*(uy+(T)1));

    return invRho*invRho*jSqr;
}



 


static T ‪invRho(T ‪rhoBar) {
         return (T)1 / (‪rhoBar + (T)1);
     }

template<typename T, template<typename U> class Descriptor>
void BGKdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics& statistics )
{
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T uSqr = dynamicsTemplates<T,Descriptor>::guo_bgk_collision(cell, rhoBar, j, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

static T ‪guo_bgk_collision(‪Array& f, T rhoBar, ‪Array const& j, T omega) 
{
     T invRho = Descriptor::invRho(rhoBar);
     const T jSqr = ‪VectorTemplateImpl::normSqr(j);
     for (‪plint iPop=0; iPop < Descriptor::q; ++iPop) {
         f[iPop] *= (T)1-omega;
         f[iPop] += omega * ‪dynamicsTemplatesImpl::bgk_ma2_equilibrium (
                                 iPop, rhoBar, invRho, j, jSqr );
     }
     return jSqr*invRho*invRho;
 }

