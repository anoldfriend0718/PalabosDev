#pragma once


#include "basicDynamics/isoThermalDynamics.h"
#include "basicDynamics/isoThermalDynamics.hh"
#include "basicDynamics/comprehensiveIsoThermalDynamics.h"
#include "core/dynamicsIdentifiers.h"
#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "plog/Log.h"
#include <map>
#include <memory>
#include <plog/Severity.h>
#include <string>

namespace plb {
enum DynamicsName { BGK_Ma2 = 0, RM, HM, CM, CHM, K, GH, RR };

enum HOOmega { SRT = 0, REG };

template <typename T, template <typename U> class Descriptor>
class LBMModelParser2D {
public:
  LBMModelParser2D(DynamicsName dynamics, HOOmega hoOmega, T omega_);
  LBMModelParser2D(std::string dynamicsName, std::string hoOmegaName, T omega_);
  LBMModelParser2D() = default;

  std::unique_ptr<Dynamics<T, Descriptor>> getDynamics() const {
    std::unique_ptr<Dynamics<T, Descriptor>> dynamics(
        isoThermalDynamics->clone());
    return dynamics;
  }

  Array<T, Descriptor<T>::numRelaxationTimes> getAllOmega() const {
    return allOmega;
  }

private:
  void setDynamics(DynamicsName &dynamics);
  void setOmega(HOOmega hoOmega);
  void setAllHOOmega();

  static const std::map<std::string, DynamicsName> dynamicsNameMap;
  static const std::map<std::string, HOOmega> hoOmegaMap;

  std::unique_ptr<IsoThermalBulkDynamics<T, Descriptor>> isoThermalDynamics;
  DynamicsName dynamicsName;
  T omega;
  Array<T, Descriptor<T>::numRelaxationTimes> allOmega;
};

} // namespace plb
