#ifndef POROSITY_EXTERNAL_FORCE_TEMPLATES_2D_H
#define POROSITY_EXTERNAL_FORCE_TEMPLATES_2D_H

#include "core/globalDefs.h"

namespace plb {
    
template<typename T>
struct Porosity_externalForceTemplatesImpl<T, descriptors::ForcedD2Q9Descriptor<T> > 
{
static void Porosity_addGuoForce (
                Array<T,descriptors::ForcedD2Q9Descriptor<T>::q>& f,
                Array<T,descriptors::ForcedD2Q9Descriptor<T>::d> const& u, T omega, T amplitude ,T porosity, T KVC, T k0 )
{
    T force[2];
    
    T invK = (((T)1-porosity)*((T)1-porosity)) / (k0 * (porosity*porosity*porosity));
    
    //T Fp = (T)1.75 / (pow(150,0.5)*pow(porosity,1.5));

    force[0] = - porosity * KVC * u[0] * invK ;//- porosity * Fp * u[0]* fabs(u[0]) / (pow(K,0.5));//Guo不带孔隙率,括号中是(c-u)/cs^2+(cu)c/cs^4
    force[1] = - porosity * KVC * u[1] * invK ;//- porosity * Fp * u[1]* fabs(u[1]) / (pow(K,0.5));

    T mu = amplitude*((T)1-omega/(T)2);
    
    static const T oneOver3 = (T)1/(T)3;
    static const T oneOver12 = (T)1/(T)12;
    static const T fourOver3 = (T)4/(T)3;
    
    const T twoUx   = (T)2*u[0];
    const T threeUx = (T)3*u[0];
    
    const T twoUy   = (T)2*u[1];
    const T threeUy = (T)3*u[1];

    f[0] += -fourOver3*mu*(force[0]*u[0]+force[1]*u[1])/porosity;
    
    f[1] += oneOver12*mu*(force[0]*(-porosity+twoUx-threeUy)+force[1]*(porosity+twoUy-threeUx))/porosity;
    
    f[2] += oneOver3*mu*(force[0]*(-porosity+twoUx)-force[1]*u[1])/porosity;
    
    f[3] += oneOver12*mu*(force[0]*(-porosity+twoUx+threeUy)+force[1]*(-porosity+twoUy+threeUx))/porosity;
    
    f[4] += -oneOver3*mu*(force[0]*u[0]+force[1]*(porosity-twoUy))/porosity;
    
    f[5] += oneOver12*mu*(force[0]*(porosity+twoUx-threeUy)+force[1]*(-porosity+twoUy-threeUx))/porosity;
    
    f[6] += oneOver3*mu*(force[0]*(porosity+twoUx)-force[1]*u[1])/porosity;
    
    f[7] += oneOver12*mu*(force[0]*(porosity+twoUx+threeUy)+force[1]*(porosity+twoUy+threeUx))/porosity;
    
    f[8] += -oneOver3*mu*(force[0]*u[0]+force[1]*(-porosity-twoUy))/porosity;
}

};

}  // namespace plb

#endif
