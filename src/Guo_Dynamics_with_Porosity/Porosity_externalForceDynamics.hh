#ifndef POROSITY_EXTERNAL_FORCE_DYNAMICS_HH
#define POROSITY_EXTERNAL_FORCE_DYNAMICS_HH

#include "basicDynamics/isoThermalDynamics.h"
#include "complexDynamics/smagorinskyDynamics.hh"
#include "core/cell.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/offEquilibriumTemplates.h"
#include "latticeBoltzmann/d3q13Templates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/externalForceTemplates.h"
#include "core/latticeStatistics.h"
#include <algorithm>
#include "core/dynamicsIdentifiers.h"
#include <limits>

#include "basicDynamics/externalForceDynamics.h"
#include "basicDynamics/externalForceDynamics.hh"
#include "Porosity_externalForceDynamics.h"
#include "Porosity_addDynamicParams.h"
#include "Porosity_externalForceTemplates.h"

namespace plb {


/* *************** Class Porosity_GuoExternalForceBGKdynamics ********************************** */

template<typename T, template<typename U> class Descriptor>
T Porosity_GuoExternalForceBGKdynamics<T,Descriptor>::getParameter(plint whichParameter) const {
    switch (whichParameter) {
        case dynamicParams::porosity   : return this->getPorosity();
        //case dynamicParams::KVC        : return this->getKVC();
        //case dynamicParams::k0         : return this->getk0();
    };
    return 0.;
}

template<typename T, template<typename U> class Descriptor>
void Porosity_GuoExternalForceBGKdynamics<T,Descriptor>::setParameter(plint whichParameter, T value) {
    switch (whichParameter) {
        case dynamicParams::porosity   : setPorosity(value);
        //case dynamicParams::KVC        : setKVC(value);
        //case dynamicParams::k0         : setk0(value);
    };
}

template<typename T, template<typename U> class Descriptor>
T Porosity_GuoExternalForceBGKdynamics<T,Descriptor>::getPorosity() const {
    return porosity;
}

template<typename T, template<typename U> class Descriptor>
void Porosity_GuoExternalForceBGKdynamics<T,Descriptor>::setPorosity(T porosity_) {
    porosity = porosity_;
}

template<typename T, template<typename U> class Descriptor>
T Porosity_GuoExternalForceBGKdynamics<T,Descriptor>::getKVC() const {
    return KVC;
}

template<typename T, template<typename U> class Descriptor>
void Porosity_GuoExternalForceBGKdynamics<T,Descriptor>::setKVC(T KVC_) {
    KVC = KVC_;
}

template<typename T, template<typename U> class Descriptor>
T Porosity_GuoExternalForceBGKdynamics<T,Descriptor>::getk0() const {
    return k0;
}

template<typename T, template<typename U> class Descriptor>
void Porosity_GuoExternalForceBGKdynamics<T,Descriptor>::setk0(T k0_) {
    k0 = k0_;
}

template<typename T, template<typename U> class Descriptor>
int Porosity_GuoExternalForceBGKdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,Porosity_GuoExternalForceBGKdynamics<T,Descriptor> >("Porosity_BGK_ExternalForce_Guo");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
Porosity_GuoExternalForceBGKdynamics<T,Descriptor>::Porosity_GuoExternalForceBGKdynamics(T omega_, T porosity_, T KVC_, T k0_ )
    : ExternalForceDynamics<T,Descriptor>(omega_),porosity(porosity_),KVC(KVC_),k0(k0_)
{ }

template<typename T, template<typename U> class Descriptor>
Porosity_GuoExternalForceBGKdynamics<T,Descriptor>::Porosity_GuoExternalForceBGKdynamics(HierarchicUnserializer& unserializer)
    : ExternalForceDynamics<T,Descriptor>(T())
{
    this->unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
Porosity_GuoExternalForceBGKdynamics<T,Descriptor>* Porosity_GuoExternalForceBGKdynamics<T,Descriptor>::clone() const {
    return new Porosity_GuoExternalForceBGKdynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int Porosity_GuoExternalForceBGKdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void Porosity_GuoExternalForceBGKdynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    Porosity_GuoExternalForceBGKdynamics<T,Descriptor>::serialize(serializer);
    serializer.addValue(porosity);
    serializer.addValue(KVC);
    serializer.addValue(k0);
}

template<typename T, template<typename U> class Descriptor>
void Porosity_GuoExternalForceBGKdynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    Porosity_GuoExternalForceBGKdynamics<T,Descriptor>::unserialize(unserializer);
    unserializer.readValue(porosity);
    unserializer.readValue(KVC);
    unserializer.readValue(k0);
}

template<typename T, template<typename U> class Descriptor>
void Porosity_GuoExternalForceBGKdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics& statistics )
{
    //calculate some parms
    T porosity = this->getPorosity();
    /*T invK =  (((T)1-porosity)*((T)1-porosity)) / k0 * (porosity*porosity*porosity);
    T Fp = (T)1.75 / (pow(150,0.5)*pow(porosity,1.5));
    
    T c0 = ((T)1 + porosity * KVC * invK / (T)2 ) / (T)2;
    T c1 = (porosity * Fp * pow(invK,0.5)) / (T)2;*/
    
    T rhoBar = this->computeRhoBar(cell);//计算f的总和=Rho-1
    Array<T,Descriptor<T>::d> u, j;
    this->computeVelocity(cell, u);//计算速度

    T rho = Descriptor<T>::fullRho(rhoBar);//计算密度
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
    {
        j[iD] = rho * u[iD];//此时j=Rho*u
    }
    T uSqr = Porosity_dynamicsTemplates<T,Descriptor>::Porosity_bgk_collision(cell, rhoBar, j, this->getOmega(), porosity);
    //Porosity_externalForceTemplates<T,Descriptor>::Porosity_addGuoForce(cell, u, this->getOmega(), (T)1, porosity, KVC, k0);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}
    
template<typename T, template<typename U> class Descriptor>
void Porosity_GuoExternalForceBGKdynamics<T,Descriptor>::collideExternal (
        Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat)
{
    //calculate some parms
    T porosity = this->getPorosity();
    /*T invK =  (((T)1-porosity)*((T)1-porosity)) / k0 * (porosity*porosity*porosity);
    T Fp = (T)1.75 / (pow(150,0.5)*pow(porosity,1.5));
    
    T c0 = ((T)1 + porosity * KVC * invK / (T)2 ) / (T)2;
    T c1 = (porosity * Fp * pow(invK,0.5)) / (T)2;*/
    
    Array<T,Descriptor<T>::d> u;
    this->computeVelocityExternal(cell, rhoBar, j, u);
    T rho = Descriptor<T>::fullRho(rhoBar);
    Array<T,Descriptor<T>::d> newJ;
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
    {
        newJ[iD] = rho * u[iD];
    }
    
    T uSqr = Porosity_dynamicsTemplates<T,Descriptor>::Porosity_bgk_collision(cell, rhoBar, newJ, this->getOmega(), porosity);
   
    //Porosity_externalForceTemplates<T,Descriptor>::Porosity_addGuoForce(cell, u, this->getOmega(), (T)1, porosity, KVC, k0);

    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
T Porosity_GuoExternalForceBGKdynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return Porosity_dynamicsTemplates<T,Descriptor>::Porosity_bgk_equilibrium(iPop, rhoBar, invRho, j, jSqr, this->getPorosity());
}

//compute Velocity in this dynamic (not use pubic member functions inherited from plb::ExternalForceDynamics ) 
template<typename T, template<typename U> class Descriptor>
void Porosity_GuoExternalForceBGKdynamics<T,Descriptor>::computeVelocity (
        Cell<T,Descriptor> const& cell, Array<T,Descriptor<T>::d>& u) const
{
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell,rhoBar, j);
    
    T porosity= this->getPorosity();
    T KVC= this->getKVC();
    T k0= this->getk0();
    T invK = (((T)1-porosity)*((T)1-porosity)) / (k0 * (porosity*porosity*porosity));
    T Fp = (T)1.75 / (pow(150,0.5)*pow(porosity,1.5));
    T invRho = Descriptor<T>::invRho(rhoBar);

    T c0 = ((T)1 + porosity * KVC * invK / (T)2 ) / (T)2;
    T c1 = (porosity * Fp * pow(invK,0.5)) / (T)2;
    
    T v[2];
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
    {
        v[iD] = j[iD]*invRho;
    }
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
    {
        u[iD] = v[iD] / (c0 + pow(c0*c0+c1*pow(v[0]*v[0]+v[1]*v[1],0.5),0.5));
    }
}

template<typename T, template<typename U> class Descriptor>
void Porosity_GuoExternalForceBGKdynamics<T,Descriptor>::computeVelocityExternal (
        Cell<T,Descriptor> const& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        Array<T,Descriptor<T>::d>& u ) const
{
    T porosity= this->getPorosity();
    T KVC= this->getKVC();
    T k0= this->getk0();
    T invK = (((T)1-porosity)*((T)1-porosity)) / (k0 * (porosity*porosity*porosity));
    T Fp = (T)1.75 / (pow(150,0.5)*pow(porosity,1.5));
    T invRho = Descriptor<T>::invRho(rhoBar);
    
    T c0 = ((T)1 + porosity * KVC * invK / (T)2 ) / (T)2;
    T c1 = (porosity * Fp * pow(invK,0.5)) / (T)2;

    T v[2];
    
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
    {
        v[iD] = j[iD]*invRho;
    }
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
    {
        u[iD] = v[iD] / (c0 + pow(c0*c0+c1*pow(v[0]*v[0]+v[1]*v[1],0.5),0.5));
    }
}


}  // namespace plb

#endif  // POROSITY_EXTERNAL_FORCE_DYNAMICS_HH
