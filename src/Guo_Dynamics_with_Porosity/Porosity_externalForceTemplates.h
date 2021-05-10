#ifndef POROSITY_EXTERNAL_FORCE_TEMPLATES_H
#define POROSITY_EXTERNAL_FORCE_TEMPLATES_H

#include "core/globalDefs.h"
#include "core/cell.h"
#include "core/util.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

template<typename T, class Descriptor> struct Porosity_externalForceTemplatesImpl;

template<typename T, template<typename U> class Descriptor>
struct Porosity_externalForceTemplates {

/// Add a force term after BGK collision, according to the Guo algorithm
//  make the assumption that F = -porosity*v*u/K, K = k0*porosity^3/(1-porosity)^2
static void Porosity_addGuoForce( Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::d> const& u,
                         T omega, T amplitude, T porosity, T KVC, T k0)
{
    Porosity_externalForceTemplatesImpl<T, Descriptor<T> >
        ::Porosity_addGuoForce(cell.getRawPopulations(), u, omega, amplitude, porosity, KVC, k0);
}

};

template<typename T, class Descriptor>
struct Porosity_externalForceTemplatesImpl {

static void Porosity_addGuoForce( Array<T,Descriptor::q>& f,
                         Array<T,Descriptor::d> const& u,
                         T omega, T amplitude, T porosity, T KVC, T k0)
{
    T force[2];
    T invK = (((T)1-porosity)*((T)1-porosity)) / (k0 * (porosity*porosity*porosity));
    //T Fp = (T)1.75 / (pow(150,0.5)*pow(porosity,1.5));
    for (int iD=0; iD < Descriptor::d; ++iD) 
    {
        force[iD] =  - porosity * KVC * u[iD] * invK ;//- porosity * Fp * u[iD]* fabs(u[iD]) / (pow(K,0.5));//Guo不带孔隙率,括号中是(c-u)/cs^2+(cu)c/cs^4
    }     

    for (plint iPop=0; iPop < Descriptor::q; ++iPop) 
    {
        T c_u = T();
        for (int iD=0; iD<Descriptor::d; ++iD) 
        {
            c_u += Descriptor::c[iPop][iD]*u[iD];
        }
        c_u *= Descriptor::invCs2 * Descriptor::invCs2;
        T forceTerm = T();
        for (int iD=0; iD < Descriptor::d; ++iD) 
        {
            forceTerm +=
                (   ((T)Descriptor::c[iPop][iD]-u[iD]/porosity) * Descriptor::invCs2
                     + c_u * (T)Descriptor::c[iPop][iD]/porosity
                )
                * force[iD];//Guo不带孔隙率,括号中是(c-u)/cs^2+(cu)c/cs^4
        }                   //Guo带孔隙率,括号中是(c-u/porosity)/cs^2+(cu)c/(cs^4*porosity)
        forceTerm *= Descriptor::t[iPop];
        forceTerm *= 1-omega/(T)2;
        forceTerm *= amplitude;
        f[iPop] += forceTerm;
    }
}

};  // struct externalForceTemplates

}  // namespace plb

#include "Porosity_externalForceTemplates2D.h"

#endif
