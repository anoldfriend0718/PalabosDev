#ifndef RESIDUAL_TRACER_HH
#define RESIDUAL_TRACER_HH

#include "core/globalDefs.h"
#include "core/util.h"
#include "customizedutil/residualTracer.h"
#include "dataProcessors/dataAnalysisWrapper2D.h"
#include "dataProcessors/dataAnalysisWrapper2D.hh"
#include "io/parallelIO.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>

namespace plb {

namespace util {

/////////// Class ResidualTracer2D ////////////////////////

template <typename T>
ResidualTracer2D<T>::ResidualTracer2D(T u, T L, plint _nx, plint _ny,
                                      T _epsilon)
    : deltaT((plint)(L / u / 2.)), epsilon(_epsilon), nx(_nx), ny(_ny), t(0),
      converged(false), absoluteResidual(MultiScalarField2D<T>(_nx, _ny, 0.0)) {
}

template <typename T> plint ResidualTracer2D<T>::getDeltaT() const {
  return deltaT;
}

template <typename T>
void ResidualTracer2D<T>::measure(MultiScalarField2D<T> &currentField,
                                  MultiScalarField2D<T> &previousField,
                                  Box2D &domain, bool doPrint) {
  subtract(currentField, previousField, absoluteResidual, domain);
  T relativeError = computeSum(*computePower(absoluteResidual, 2.0, domain)) /
                    computeSum(*computePower(currentField, 2.0, domain));

  values.push_back(relativeError);
  if ((plint)values.size() > std::abs(deltaT)) {
    values.erase(values.begin());
    if (doPrint && t % deltaT == 0) {
      T average = computeAverage();
      pcout << "average relative error=" << average << std::endl;
    }
  }
  ++t;
}

template <typename T> void ResidualTracer2D<T>::resetScale(T u, T L) {
  t = t % deltaT;
  deltaT = (plint)(L / u / 2.);
  if ((plint)values.size() > std::abs(deltaT)) {
    values.erase(values.begin(), values.begin() + (values.size() - deltaT));
  }
}

template <typename T> void ResidualTracer2D<T>::resetValues() {
  t = 0;
  if ((plint)values.size() > 0) {
    values.erase(values.begin(), values.begin() + values.size());
  }
}

template <typename T> bool ResidualTracer2D<T>::hasConverged() const {
  if ((plint)values.size() < std::abs(deltaT)) {
    return false;
  } else {
    T average = computeAverage();
    if (!util::isNaN(average)) {
      bool isConvergence = average < epsilon;
      if (!isConvergence) {
        return false;
      }
      pcout << "simulation is converged with the average residual error: "
            << average << std::endl;
      return true;
    }

    else {
      pcout << "simulation diverged.\n";
      return true;
    }
  }
}

template <typename T> T ResidualTracer2D<T>::computeAverage() const {
  return accumulate(values.begin(), values.end(), 0.) / values.size();
}

template <typename T> void ResidualTracer2D<T>::setEpsilon(T epsilon_) {
  epsilon = epsilon_;
}

} // namespace util

} // namespace plb

#endif
