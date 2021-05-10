#ifndef MT_DYNAMICS_TEMPLATES_2D_H
#define MT_DYNAMICS_TEMPLATES_2D_H

namespace plb {

/// All helper functions are inside this structure
template<typename T>
struct MTdynamicsTemplatesImpl<T,descriptors::D2Q5DescriptorBase<T> >
{
    
typedef descriptors::D2Q5DescriptorBase<T> Descriptor;

static T MTequilibrium(plint iPop, T rhoBar, Array<T,Descriptor::d> const& jEq, T porosity) 
{
    return Descriptor::t[iPop] * (rhoBar + Descriptor::invCs2 * 
            (Descriptor::c[iPop][0]*jEq[0]+Descriptor::c[iPop][1]*jEq[1]) / porosity);
}


static T MTcollision (
        Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& jEq, 
        T omega, T porosity ) 
{
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = jEq[0]*jEq[0] + jEq[1]*jEq[1];
    
    const T oneMinusOmega = (T)1 - omega;
    const T halfOmega = (T)0.5 * omega;
    const T cs2RhoBar = Descriptor::cs2*rhoBar;
    
    f[0] = oneMinusOmega*f[0]+omega*((T)1-2*Descriptor::cs2)*rhoBar;
    
    f[1] = oneMinusOmega*f[1]+halfOmega*(cs2RhoBar-jEq[0]/porosity);
    f[2] = oneMinusOmega*f[2]+halfOmega*(cs2RhoBar-jEq[1]/porosity);
    f[3] = oneMinusOmega*f[3]+halfOmega*(cs2RhoBar+jEq[0]/porosity);
    f[4] = oneMinusOmega*f[4]+halfOmega*(cs2RhoBar+jEq[1]/porosity);
    
    return jSqr*invRho*invRho;
}

};  // struct advectionDiffusionDynamicsTemplatesImpl

template<typename T>
struct MTdynamicsTemplatesImpl<T,descriptors::D2Q9DescriptorBase<T> >
{
    
typedef descriptors::D2Q9DescriptorBase<T> D;

static T MTequilibrium(plint iPop, T rhoBar, Array<T,D::d> const& jEq, T porosity) 
{
    return Porosity_dynamicsTemplatesImpl<T,D>::Porosity_bgk_equilibrium(iPop, rhoBar, D::invRho(rhoBar), jEq, VectorTemplateImpl<T,D::d>::normSqr(jEq), porosity );
}

static T MTcollision ( Array<T,D::q>& f, T rhoBar, Array<T,D::d> const& jEq, 
        T omega, T porosity ) 
{
    T invRho = D::invRho(rhoBar);
    T jSqr = Porosity_dynamicsTemplatesImpl<T,D>::Porosity_bgk_collision(f, rhoBar, jEq, omega, porosity);
    
    return jSqr*invRho*invRho;
}
  
};  // struct advectionDiffusionDynamicsTemplatesImpl

}  // namespace plb

#endif
