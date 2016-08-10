// -*- C++ -*-
/**
 * \file test_unsaturated_zone_tracer_solver.cc
 * \author Jeffrey Rutledge jrutledge@hmc.edu
 *
 * \brief Tests for correctness of methods in unsaturated_zone_tracer_solver.h
 */

#include "unsaturated_zone_tracer_solver.h"

#include <array>

#include <gtest/gtest.h>

namespace solver = unsaturated_zone_tracer_solver;
namespace internal = unsaturated_zone_tracer_solver_internal;

//--------------------------------------------------
//           TESTS
//--------------------------------------------------

TEST(ThomasAlgorithim, SimpleCase) {
  constexpr size_t MATRIX_SIZE = 4;
  std::array<double, MATRIX_SIZE> rhs_vector = {1., 0., 0., 1.};
  std::array<double, MATRIX_SIZE> solution = internal::ThomasAlgorithimSingleValue<MATRIX_SIZE>(
      -1., 2., -1., rhs_vector);

  std::array<double, MATRIX_SIZE> true_solution = {1., 1., 1., 1.};
  for (size_t i = 0; i < MATRIX_SIZE; ++i) {
    ASSERT_DOUBLE_EQ(true_solution[i], solution[i]);
  }
}

//--------------------------------------------------
//           RUNNING THE TESTS
//--------------------------------------------------

// Called if the test runs too long.
static void timeout_handler(int) {
  // We go super-low-level here, because we can't trust anything in
  // the C/C++ library to really be working right.
  write(STDERR_FILENO, "Timeout occurred!\n", 18);
  abort();
}

/// Run tests
int main(int argc, char** argv) {
  // Initalize testing environment
  ::testing::InitGoogleTest(&argc, argv);

  const size_t MAX_RUNTIME = 30;  // Timeout, in seconds
  if (not::testing::GTEST_FLAG(break_on_failure)) {
    signal(SIGALRM, timeout_handler);  // What to call when timer expires
    alarm(MAX_RUNTIME);                // set timer at MAX_RUNTIME seconds
  }
  return RUN_ALL_TESTS();
}
