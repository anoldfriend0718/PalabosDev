#pragma once

#include "LBMModelParser2D.h"

namespace plb {
template <typename T, template <typename U> class Descriptor>
LBMModelParser2D<T, Descriptor>::LBMModelParser2D(DynamicsName dynamics,
                                                  HOOmega hoOmega, T omega_)
    : dynamicsName(dynamics), omega(omega_) {
  setDynamics(dynamics);
  setOmega(hoOmega);
  setAllHOOmega();
};

template <typename T, template <typename U> class Descriptor>
LBMModelParser2D<T, Descriptor>::LBMModelParser2D(std::string dynamicsName,
                                                  std::string hoOmegaName,
                                                  T omega_)
    : LBMModelParser2D(dynamicsNameMap.at(dynamicsName),
                       hoOmegaMap.at(hoOmegaName), omega_){};

template <typename T, template <typename U> class Descriptor>
void LBMModelParser2D<T, Descriptor>::setDynamics(DynamicsName &dynamics) {
  switch (dynamics) {
  case (BGK_Ma2):
    isoThermalDynamics.reset(new BGKdynamics<T, Descriptor>(omega));
    break;
  case (RM):
    isoThermalDynamics.reset(new RMdynamics<T, Descriptor>(omega));
    break;
  case (HM):
    isoThermalDynamics.reset(new HMdynamics<T, Descriptor>(omega));
    break;
  case (CM):
    isoThermalDynamics.reset(new CMdynamics<T, Descriptor>(omega));
    break;
  case (CHM):
    isoThermalDynamics.reset(new CHMdynamics<T, Descriptor>(omega));
    break;
  case (GH):
    isoThermalDynamics.reset(new GHdynamics<T, Descriptor>(omega));
    break;
  case (RR):
    isoThermalDynamics.reset(new RRdynamics<T, Descriptor>(omega));
    break;
  case (K):
    isoThermalDynamics.reset(new Kdynamics<T, Descriptor>(omega));
    break;
  default:
    PLOG(plog::error) << "Error: dynamics name does not exist." << std::endl;
    exit(-1);
    break;
  }
}

template <typename T, template <typename U> class Descriptor>
void LBMModelParser2D<T, Descriptor>::setOmega(HOOmega hoOmega) {
  PLB_PRECONDITION(Descriptor<T>::numRelaxationTimes == 4);
  switch (hoOmega) {
  case (SRT):
    allOmega[0] = omega;
    allOmega[1] = omega;
    allOmega[2] = omega;
    allOmega[3] = omega;
    break;
  case (REG):
    allOmega[0] = omega;
    allOmega[1] = omega;
    allOmega[2] = 1;
    allOmega[3] = 1;
    break;
  default:
    PLOG(plog::error) << "Error: high order omega catelog does not exist."
                      << std::endl;
    exit(-1);
    break;
  }
}

template <typename T, template <typename U> class Descriptor>
void LBMModelParser2D<T, Descriptor>::setAllHOOmega() {
  switch (dynamicsName) {
  case (BGK_Ma2):
    // not high order omega for BGK_Ma2
    break;
  case (RM):
    RMdynamics<T, Descriptor>::allOmega = allOmega;
    break;
  case (HM):
    HMdynamics<T, Descriptor>::allOmega = allOmega;
    break;
  case (CM):
    CMdynamics<T, Descriptor>::allOmega = allOmega;
    break;
  case (CHM):
    CHMdynamics<T, Descriptor>::allOmega = allOmega;
    break;
  case (GH):
    GHdynamics<T, Descriptor>::allOmega = allOmega;
    break;
  case (RR):
    RRdynamics<T, Descriptor>::allOmega = allOmega;
    break;
  case (K):
    Kdynamics<T, Descriptor>::allOmega = allOmega;
    break;
  default:
    PLOG(plog::error) << "Error: dynamics name does not exist." << std::endl;
    exit(-1);
    break;
  }
}

template <typename T, template <typename U> class Descriptor>
const std::map<std::string, DynamicsName>
    LBMModelParser2D<T, Descriptor>::dynamicsNameMap = {
        {"BGK_Ma2", BGK_Ma2}, {"RM", RM}, {"HM", HM}, {"CM", CM},
        {"CHM", CHM},         {"K", K},   {"GH", GH}, {"RR", RR}};

template <typename T, template <typename U> class Descriptor>
const std::map<std::string, HOOmega>
    LBMModelParser2D<T, Descriptor>::hoOmegaMap = {{"SRT", SRT}, {"REG", REG}};

} // namespace plb
