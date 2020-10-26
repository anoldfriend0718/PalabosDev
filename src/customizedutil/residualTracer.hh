#ifndef RESIDUAL_TRACER_HH
#define RESIDUAL_TRACER_HH

#include "core/globalDefs.h"
#include "customizedutil/residualTracer.h"
#include "dataProcessors/dataAnalysisWrapper2D.h"
#include "dataProcessors/dataAnalysisWrapper2D.hh"
#include "io/parallelIO.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>

namespace plb {

namespace util {

/////////// Class ResidualTracer2D ////////////////////////

template <typename T>
ResidualTracer2D<T>::ResidualTracer2D(plint _count, plint _nx, plint _ny,
                                      T _epsilon)
    : count(std::abs(_count)), epsilon(_epsilon), nx(_nx), ny(_ny), t(0),
      converged(false),
      absoluteResidualField(MultiScalarField2D<T>(_nx, _ny, 0.0)) {}

template <typename T>
void ResidualTracer2D<T>::measure(MultiScalarField2D<T> &currentField,
                                  MultiScalarField2D<T> &previousField,
                                  Box2D domain, bool doPrint) {
  subtract(currentField, previousField, absoluteResidualField, domain);
  T relativeError =
      sqrt(computeSum(*computePower(absoluteResidualField, 2.0, domain)) /
           computeSum(*computePower(currentField, 2.0, domain)));

  relativeErrors.push_back(relativeError);
  if ((plint)relativeErrors.size() > count) {
    relativeErrors.erase(relativeErrors.begin());
    if (doPrint && t % count == 0) {
      T average = computeAverage();
      pcout << "average relative error=" << average << std::endl;
    }
  }
  ++t;
}

template <typename T> void ResidualTracer2D<T>::resetCount(plint _count) {
  t = t % count;
  count = std ::abs(_count);
  if ((plint)relativeErrors.size() > count) {
    relativeErrors.erase(relativeErrors.begin(),
                         relativeErrors.begin() +
                             (relativeErrors.size() - count));
  }
}

template <typename T> void ResidualTracer2D<T>::resetValues() {
  t = 0;
  if ((plint)relativeErrors.size() > 0) {
    relativeErrors.erase(relativeErrors.begin(),
                         relativeErrors.begin() + relativeErrors.size());
  }
}

template <typename T> bool ResidualTracer2D<T>::hasConverged() const {
  if ((plint)relativeErrors.size() < count) {
    return false;
  }

  T average = computeAverage();
  if (util::isNaN(average)) {
    pcout << "simulation diverged.\n";
    return true;
  }

  bool isConvergence = average < epsilon;
  if (!isConvergence) {
    return false;
  }
  pcout << std::endl
        << "simulation is converged with the average residual error: "
        << average << std::endl;
  return true;
} // namespace util

template <typename T> T ResidualTracer2D<T>::computeAverage() const {

  return accumulate(relativeErrors.begin(), relativeErrors.end(), 0.) /
         relativeErrors.size();
}

template <typename T> void ResidualTracer2D<T>::setEpsilon(T epsilon_) {
  epsilon = epsilon_;
}

} // namespace util

} // namespace plb

#endif
