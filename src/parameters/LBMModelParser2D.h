#include "basicDynamics/comprehensiveIsoThermalDynamics.h"
#include "basicDynamics/isoThermalDynamics.h"
#include "core/dynamicsIdentifiers.h"
#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "plog/Log.h"
#include <memory>
#include <plog/Severity.h>

namespace plb {
enum DynamicsName { BGK_Ma2 = 0, RM, HM, CM, CHM, K, GH, RR };

enum HOOmega { SRT = 0, REG };

template <typename T, template <typename U> class Descriptor>
class LBMModelParser2D {
public:
  LBMModelParser2D(DynamicsName dynamics, HOOmega hoOmega, T omega_)
      : dynamicsName(dynamics),omega(omega_) {
    setDynamics(dynamics);
  setOmega(hoOmega);
  }

  LBMModelParser2D() = default;

  std::unique_ptr<Dynamics<T, Descriptor>> getDynamics() const {
    std::unique_ptr<Dynamics<T, Descriptor>> dynamics(
        isoThermalDynamics->clone());
    return dynamics;
  }

  Array<T, Descriptor<T>::numRelaxationTimes> getAllOmega() const {
    return allOmega;
  }

  void setAllOmega()
  {
    switch (dynamicsName) {
    case (BGK_Ma2):
      BGKdynamics<T, Descriptor>::allOmega = allOmega;
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

private:
  std::unique_ptr<IsoThermalBulkDynamics<T, Descriptor>> isoThermalDynamics;
  DynamicsName dynamicsName;
  T omega;
  Array<T, Descriptor<T>::numRelaxationTimes> allOmega;

  void setDynamics(DynamicsName &dynamics) {
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
  void setOmega(HOOmega hoOmega) {
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
};
} // namespace plb
