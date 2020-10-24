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
  /** \param count relative error measurement count for
   *  averaging the output relative error
   * \param _epsilon Precision of the convergence.
   * \param nx domain Nx
   * \param ny domain Ny
   */
  ResidualTracer2D(plint _count, plint _nx, plint _ny, T epsilon);
  /// Change count
  void resetCount(plint _count);
  /// reinitializes the values
  void resetValues();
  /// Feed the object with a new measured scalar.
  void measure(MultiScalarField2D<T> &currentField,
               MultiScalarField2D<T> &previousField, Box2D &domain,
               bool doPrint = false);
  /// Test for convergence, with respect to stdDev.
  bool hasConverged() const;

  void setEpsilon(T epsilon_);

private:
  plint count;
  T epsilon;
  plint nx;
  plint ny;
  plint t;
  bool converged;
  std::deque<T> relativeErrors;
  MultiScalarField2D<T> absoluteResidualField;

  T computeAverage() const;
};
} // namespace util

} // namespace plb

#endif
