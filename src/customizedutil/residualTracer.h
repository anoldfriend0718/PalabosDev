#ifndef RESIDUAL_TRACER_H
#define RESIDUAL_TRACER_H

#include "core/globalDefs.h"
#include "multiBlock/multiDataField2D.h"
#include <deque>

namespace plb {

namespace util {

/// Check time-convergence of a scalar.
/** This class is useful, for example to check convergence of
 * the velocity field for the simulation of a stationary flow.
 * Convergence is claimed when the relative residual error is smaller than
 * epsilon . The statistics are taken over a macroscopic time scale of the
 * system.
 */
template <typename T> class ResidualTracer2D {
public:
  /// The only constructor.
  /** \param u The characteristic velocity of the system, for
   *          computation of the characteristic time scale.
   * \param L The characteristic length of the system, for
   *          computation of the characteristic time scale.
   * \param _epsilon Precision of the convergence.
   * \param nx domain Nx
   * \param ny domain Ny
   */
  ResidualTracer2D(T u, T L, plint _nx, plint _ny, T epsilon);
  /// Change values of u and L to update characteristic scales of the system.
  void resetScale(T u, T L);
  /// reinitializes the values
  void resetValues();
  /// Get characteristic time scale.
  plint getDeltaT() const;
  /// Feed the object with a new measured scalar.
  void measure(MultiScalarField2D<T> &currentField,
               MultiScalarField2D<T> &previousField, Box2D &domain,
               bool doPrint = false);
  /// Test for convergence, with respect to stdDev.
  bool hasConverged() const;

  void setEpsilon(T epsilon_);

private:
  plint deltaT;
  T epsilon;
  plint nx;
  plint ny;
  plint t;
  bool converged;
  std::deque<T> values;
  MultiScalarField2D<T> absoluteResidual;

  T computeAverage() const;
};
} // namespace util

} // namespace plb

#endif
