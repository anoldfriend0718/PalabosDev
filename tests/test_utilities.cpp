#include "core/geometry2D.h"
#include "core/globalDefs.h"
#include "customizedutil/residualTracer.h"
#include "customizedutil/residualTracer.hh"
#include "dataProcessors/dataAnalysisWrapper2D.h"
#include "palabos2D.h"
#include "palabos2D.hh"
#include "gtest/gtest.h"

using namespace testing;
using namespace plb;

typedef double T;
int main(int argc, char **argv) {
  plbInit(&argc, &argv);
  InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

TEST(UtilityMethods, TestResidualTracer2D) {
  plint nx = 2;
  plint ny = 2;
  util::ResidualTracer2D<T> residualTracer(1, nx, ny, 1e-6);
  MultiScalarField2D<T> previousField(nx, ny, 1.0 - 1e-7);
  MultiScalarField2D<T> currentField(nx, ny, 1.0);
  residualTracer.measure(currentField, previousField,
                         currentField.getBoundingBox());
  auto relativeErrors = residualTracer.getRelativeErrors();
  ASSERT_EQ(relativeErrors.size(), 1);
  ASSERT_NEAR(relativeErrors[0], 1e-7, 1e-9);
  auto isConvergence = residualTracer.hasConverged();
  ASSERT_TRUE(isConvergence);
}

TEST(UtilityMethods, TestComputeSum) {
  MultiScalarField2D<T> field(2, 2, (double)1.0);
  T result = computeSum(field, field.getBoundingBox());
  ASSERT_EQ(result, 4.0);
}
