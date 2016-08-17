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

// Sample data for tests
constexpr std::array<double, 149> SAMPLE_SURFACE_CONCENTRATION = {
  0.011857622, 
  0.013833893, 
  0.015810162999999999, 
  0.018280501000000001, 
  0.020750839, 
  0.024209312, 
  0.027667786, 
  0.032114393999999997, 
  0.036561002000000002, 
  0.041995746, 
  0.047430489999999999, 
  0.057805909000000003, 
  0.068181328999999999, 
  0.088438100000000006, 
  0.108694872, 
  0.15316095599999999, 
  0.19762704, 
  0.26679650399999999, 
  0.335965968, 
  0.44466084, 
  0.553355712, 
  0.70157599199999998, 
  0.84979627199999996, 
  1.06224534, 
  1.274694408, 
  1.571134968, 
  1.8675755279999999, 
  2.2480075799999999, 
  2.6284396320000001, 
  3.0928631759999998, 
  3.55728672, 
  4.1353458119999997, 
  4.7134049039999999, 
  5.3853368399999999, 
  6.0572687759999999, 
  6.7390820639999998, 
  7.4208953519999996, 
  8.0681239080000005, 
  8.7153524640000004, 
  9.4762165679999999, 
  10.237080669999999, 
  11.23015655, 
  12.22323242, 
  13.497926830000001, 
  14.772621239999999, 
  16.348696879999999, 
  17.924772529999998, 
  19.82199211, 
  21.719211699999999, 
  23.912871840000001, 
  26.10653198, 
  28.566988630000001, 
  31.027445279999998, 
  33.78434249, 
  36.541239699999998, 
  39.668687599999998, 
  42.796135509999999, 
  46.383066290000002, 
  49.969997059999997, 
  54.080639499999997, 
  58.191281930000002, 
  62.746585199999998, 
  67.301888469999994, 
  72.356200020000003, 
  77.410511569999997, 
  83.136755050000005, 
  88.862998540000007, 
  95.251292599999999, 
  101.6395867, 
  107.7907283, 
  113.9418699, 
  121.07620609999999, 
  128.21054219999999, 
  134.7569379, 
  141.3033336, 
  148.10121079999999, 
  151.5001493, 
  154.89908790000001, 
  157.24782540000001, 
  160.80271579999999, 
  167.2725307, 
  171.62647279999999, 
  174.0226395, 
  178.4938373, 
  181.88634830000001, 
  186.85983529999999, 
  191.72627700000001, 
  196.1171085, 
  200.76353109999999, 
  202.9189283, 
  210.21809690000001, 
  214.4178938, 
  222.06355959999999, 
  226.01764639999999, 
  232.2015682, 
  238.5505043, 
  246.72283279999999, 
  252.12450079999999, 
  255.79928179999999, 
  260.19604199999998, 
  263.27104029999998, 
  262.35473300000001, 
  265.82530459999998, 
  265.93761979999999, 
  268.52844420000002, 
  270.34590789999999, 
  270.40190080000002, 
  270.86466560000002, 
  271.08106170000002, 
  270.19439770000002, 
  269.38546500000001, 
  269.49020460000003, 
  269.19838279999999, 
  267.78472679999999, 
  267.24126630000001, 
  266.45736570000003, 
  265.93959599999999, 
  264.84213519999997, 
  263.86291820000002, 
  262.14097179999999, 
  262.0365615, 
  261.1195955, 
  261.0813885, 
  259.41016519999999, 
  258.89997099999999, 
  256.99357750000001, 
  256.88488539999997, 
  254.87506959999999, 
  253.47261219999999, 
  251.99966950000001, 
  250.93415759999999, 
  249.85975260000001, 
  248.6582109, 
  247.7060022, 
  246.5212583, 
  245.19192100000001, 
  243.88498079999999, 
  242.91333929999999, 
  241.99472639999999, 
  240.6620954, 
  240.1878026, 
  238.804778, 
  237.49783780000001, 
  236.52619630000001, 
  235.60758340000001, 
  234.2749523, 
  232.8919277, 
  231.5089031, 
  230.20196290000001};

//--------------------------------------------------
//           TESTS
//--------------------------------------------------

TEST(ThomasAlgorithim, SimpleCase) {
  constexpr size_t MATRIX_SIZE = 4;
  std::array<double, MATRIX_SIZE> rhs_vector = {1., 0., 0., 1.};
  std::array<double, MATRIX_SIZE> solution =
      internal::ThomasAlgorithimSingleValue<MATRIX_SIZE>(-1., 2., -1.,
                                                         rhs_vector);

  std::array<double, MATRIX_SIZE> true_solution = {1., 1., 1., 1.};
  for (size_t i = 0; i < MATRIX_SIZE; ++i) {
    ASSERT_DOUBLE_EQ(true_solution[i], solution[i]);
  }
}

// TODO Add Tests for Fully Implicit solver
TEST(CrankNicolsonAtDepth, ExactValue) {
  // Sample parameters
  constexpr size_t TIME_STEPS = 148;
  constexpr size_t DEPTH_STEPS = 1000;
  constexpr double MAX_TIME = 74;
  constexpr double MAX_DEPTH = 200;
  constexpr double D_STAR = 9.88011475;
  constexpr double Q_STAR = 0.3;
  constexpr double THETA_STAR = 0.2255;
  constexpr double EFFECTIVE_DIFFUSION = D_STAR / THETA_STAR;
  constexpr double EFFECTIVE_VELOCITY = Q_STAR / THETA_STAR;
  constexpr double DECAY_RATE = 0.;
  constexpr double DELTA_DEPTH = MAX_DEPTH / DEPTH_STEPS;

  const auto solution = solver::CrankNicolson<TIME_STEPS, DEPTH_STEPS>(
      MAX_TIME, MAX_DEPTH, EFFECTIVE_DIFFUSION, EFFECTIVE_VELOCITY, DECAY_RATE,
      SAMPLE_SURFACE_CONCENTRATION);

  // Convert sample surface concentrations to a vector
  std::vector<double> sample_surface_concentration_vector(
      SAMPLE_SURFACE_CONCENTRATION.size());
  std::copy(SAMPLE_SURFACE_CONCENTRATION.cbegin(),
            SAMPLE_SURFACE_CONCENTRATION.cend(),
            sample_surface_concentration_vector.begin());

  // Test at depth = 0
  double requested_depth = 0;
  size_t expected_depth_step = 0;
  std::vector<double> depth_solution =
      solver::CrankNicolsonAtDepth(
          TIME_STEPS, MAX_TIME, DEPTH_STEPS, EFFECTIVE_DIFFUSION,
          EFFECTIVE_VELOCITY, DECAY_RATE, requested_depth,
          sample_surface_concentration_vector);
  // Test solution has the correct length
  ASSERT_EQ(TIME_STEPS + 1, depth_solution.size());
  // Test solution has the same values as the `CrankNicolson()` solution.
  for (size_t time_step = 0; time_step <= TIME_STEPS; ++time_step) {
    ASSERT_DOUBLE_EQ(solution[time_step][expected_depth_step],
                     depth_solution[time_step])
        << "Requested Depth: " << requested_depth << std::endl
        << "Expected Depth Step: " << expected_depth_step << std::endl
        << "Time Step: " << time_step;
  }

  // Test at middle depth
  requested_depth = 100.2;
  expected_depth_step = 100.2 / DELTA_DEPTH;
  depth_solution =
      solver::CrankNicolsonAtDepth(
          TIME_STEPS, MAX_TIME, DEPTH_STEPS, EFFECTIVE_DIFFUSION,
          EFFECTIVE_VELOCITY, DECAY_RATE, requested_depth,
          sample_surface_concentration_vector);
  // Test solution has the correct length
  ASSERT_EQ(TIME_STEPS + 1, depth_solution.size());
  // Test solution has the same values
  // Test solution has the same values as the `CrankNicolson()` solution.
  for (size_t time_step = 0; time_step <= TIME_STEPS; ++time_step) {
    ASSERT_DOUBLE_EQ(solution[time_step][expected_depth_step],
                     depth_solution[time_step])
        << "Requested Depth: " << requested_depth << std::endl
        << "Expected Depth Step: " << expected_depth_step << std::endl
        << "Time Step: " << time_step;
  }

  // Test at max depth
  requested_depth = 200;
  expected_depth_step = 200 / DELTA_DEPTH;
  depth_solution =
      solver::CrankNicolsonAtDepth(
          TIME_STEPS, MAX_TIME, DEPTH_STEPS, EFFECTIVE_DIFFUSION,
          EFFECTIVE_VELOCITY, DECAY_RATE, requested_depth,
          sample_surface_concentration_vector);
  // Test solution has the correct length
  ASSERT_EQ(TIME_STEPS + 1, depth_solution.size());
  // Test solution has the same values as the `CrankNicolson()` solution.
  for (size_t time_step = 0; time_step <= TIME_STEPS; ++time_step) {
    ASSERT_DOUBLE_EQ(solution[time_step][expected_depth_step],
                     depth_solution[time_step])
        << "Requested Depth: " << requested_depth << std::endl
        << "Expected Depth Step: " << expected_depth_step << std::endl
        << "Time Step: " << time_step;
  }
}

TEST(CrankNicolsonAtDepth, RoundDown) {
  // Sample parameters
  constexpr size_t TIME_STEPS = 148;
  constexpr size_t DEPTH_STEPS = 1000;
  constexpr double MAX_TIME = 74;
  constexpr double MAX_DEPTH = 200;
  constexpr double D_STAR = 9.88011475;
  constexpr double Q_STAR = 0.3;
  constexpr double THETA_STAR = 0.2255;
  constexpr double EFFECTIVE_DIFFUSION = D_STAR / THETA_STAR;
  constexpr double EFFECTIVE_VELOCITY = Q_STAR / THETA_STAR;
  constexpr double DECAY_RATE = 0.;
  constexpr double DELTA_DEPTH = MAX_DEPTH / DEPTH_STEPS;

  const auto solution = solver::CrankNicolson<TIME_STEPS, DEPTH_STEPS>(
      MAX_TIME, MAX_DEPTH, EFFECTIVE_DIFFUSION, EFFECTIVE_VELOCITY, DECAY_RATE,
      SAMPLE_SURFACE_CONCENTRATION);

  // Convert sample surface concentrations to a vector
  std::vector<double> sample_surface_concentration_vector(
      SAMPLE_SURFACE_CONCENTRATION.size());
  std::copy(SAMPLE_SURFACE_CONCENTRATION.cbegin(),
            SAMPLE_SURFACE_CONCENTRATION.cend(),
            sample_surface_concentration_vector.begin());

  // Test at middle depth
  double requested_depth = 100.25;
  size_t expected_depth_step = 100.2 / DELTA_DEPTH;
  std::vector<double> depth_solution =
      solver::CrankNicolsonAtDepth(
          TIME_STEPS, MAX_TIME, DEPTH_STEPS, EFFECTIVE_DIFFUSION,
          EFFECTIVE_VELOCITY, DECAY_RATE, requested_depth,
          sample_surface_concentration_vector);
  // Test solution has the correct length
  ASSERT_EQ(TIME_STEPS + 1, depth_solution.size());
  // Test solution has the same values as the `CrankNicolson()` solution.
  for (size_t time_step = 0; time_step <= TIME_STEPS; ++time_step) {
    ASSERT_DOUBLE_EQ(solution[time_step][expected_depth_step],
                     depth_solution[time_step])
        << "Requested Depth: " << requested_depth << std::endl
        << "Expected Depth Step: " << expected_depth_step << std::endl
        << "Time Step: " << time_step;
  }
}

TEST(CrankNicolsonAtDepth, RoundUp) {
  // Sample parameters
  constexpr size_t TIME_STEPS = 148;
  constexpr size_t DEPTH_STEPS = 1000;
  constexpr double MAX_TIME = 74;
  constexpr double MAX_DEPTH = 200;
  constexpr double D_STAR = 9.88011475;
  constexpr double Q_STAR = 0.3;
  constexpr double THETA_STAR = 0.2255;
  constexpr double EFFECTIVE_DIFFUSION = D_STAR / THETA_STAR;
  constexpr double EFFECTIVE_VELOCITY = Q_STAR / THETA_STAR;
  constexpr double DECAY_RATE = 0.;
  constexpr double DELTA_DEPTH = MAX_DEPTH / DEPTH_STEPS;

  const auto solution = solver::CrankNicolson<TIME_STEPS, DEPTH_STEPS>(
      MAX_TIME, MAX_DEPTH, EFFECTIVE_DIFFUSION, EFFECTIVE_VELOCITY, DECAY_RATE,
      SAMPLE_SURFACE_CONCENTRATION);

  // Convert sample surface concentrations to a vector
  std::vector<double> sample_surface_concentration_vector(
      SAMPLE_SURFACE_CONCENTRATION.size());
  std::copy(SAMPLE_SURFACE_CONCENTRATION.cbegin(),
            SAMPLE_SURFACE_CONCENTRATION.cend(),
            sample_surface_concentration_vector.begin());

  // Test at middle depth
  double requested_depth = 100.15;
  size_t expected_depth_step = 100.2 / DELTA_DEPTH;
  std::vector<double> depth_solution =
      solver::CrankNicolsonAtDepth(
          TIME_STEPS, MAX_TIME, DEPTH_STEPS, EFFECTIVE_DIFFUSION,
          EFFECTIVE_VELOCITY, DECAY_RATE, requested_depth,
          sample_surface_concentration_vector);
  // Test solution has the correct length
  ASSERT_EQ(TIME_STEPS + 1, depth_solution.size());
  // Test solution has the same values as the `CrankNicolson()` solution.
  for (size_t time_step = 0; time_step <= TIME_STEPS; ++time_step) {
    ASSERT_DOUBLE_EQ(solution[time_step][expected_depth_step],
                     depth_solution[time_step])
        << "Requested Depth: " << requested_depth << std::endl
        << "Expected Depth Step: " << expected_depth_step << std::endl
        << "Time Step: " << time_step;
  }
}

TEST(FullyImplicitAtDepth, ExactValue) {
  // Sample parameters
  constexpr size_t TIME_STEPS = 148;
  constexpr size_t DEPTH_STEPS = 1000;
  constexpr double MAX_TIME = 74;
  constexpr double MAX_DEPTH = 200;
  constexpr double D_STAR = 9.88011475;
  constexpr double Q_STAR = 0.3;
  constexpr double THETA_STAR = 0.2255;
  constexpr double EFFECTIVE_DIFFUSION = D_STAR / THETA_STAR;
  constexpr double EFFECTIVE_VELOCITY = Q_STAR / THETA_STAR;
  constexpr double DECAY_RATE = 0.;
  constexpr double DELTA_DEPTH = MAX_DEPTH / DEPTH_STEPS;

  const auto solution = solver::FullyImplicit<TIME_STEPS, DEPTH_STEPS>(
      MAX_TIME, MAX_DEPTH, EFFECTIVE_DIFFUSION, EFFECTIVE_VELOCITY, DECAY_RATE,
      SAMPLE_SURFACE_CONCENTRATION);
  //TODO Why does this print line make the test pass?
  std::cout << solution[1][501] << std::endl;

  // Convert sample surface concentrations to a vector
  std::vector<double> sample_surface_concentration_vector(
      SAMPLE_SURFACE_CONCENTRATION.size());
  std::copy(SAMPLE_SURFACE_CONCENTRATION.cbegin(),
            SAMPLE_SURFACE_CONCENTRATION.cend(),
            sample_surface_concentration_vector.begin());

  // Test at depth = 0
  double requested_depth = 0;
  size_t expected_depth_step = 0;
  std::vector<double> depth_solution =
      solver::FullyImplicitAtDepth(
          TIME_STEPS, MAX_TIME, DEPTH_STEPS, EFFECTIVE_DIFFUSION,
          EFFECTIVE_VELOCITY, DECAY_RATE, requested_depth,
          sample_surface_concentration_vector);
  // Test solution has the correct length
  ASSERT_EQ(TIME_STEPS + 1, depth_solution.size());
  // Test solution has the same values as the `FullyImplicit()` solution.
  for (size_t time_step = 0; time_step <= TIME_STEPS; ++time_step) {
    ASSERT_DOUBLE_EQ(solution[time_step][expected_depth_step],
                     depth_solution[time_step])
        << "Requested Depth: " << requested_depth << std::endl
        << "Expected Depth Step: " << expected_depth_step << std::endl
        << "Time Step: " << time_step;
  }

  // Test at middle depth
  requested_depth = 100.2;
  expected_depth_step = 100.2 / DELTA_DEPTH;
  depth_solution =
      solver::FullyImplicitAtDepth(
          TIME_STEPS, MAX_TIME, DEPTH_STEPS, EFFECTIVE_DIFFUSION,
          EFFECTIVE_VELOCITY, DECAY_RATE, requested_depth,
          sample_surface_concentration_vector);
  // Test solution has the correct length
  ASSERT_EQ(TIME_STEPS + 1, depth_solution.size());
  // Test solution has the same values as the `FullyImplicit()` solution.
  for (size_t time_step = 0; time_step <= TIME_STEPS; ++time_step) {
    ASSERT_DOUBLE_EQ(solution[time_step][expected_depth_step],
                     depth_solution[time_step])
        << "Requested Depth: " << requested_depth << std::endl
        << "Expected Depth Step: " << expected_depth_step << std::endl
        << "Time Step: " << time_step;
  }

  // Test at max depth
  requested_depth = 200;
  expected_depth_step = 200 / DELTA_DEPTH;
  depth_solution =
      solver::FullyImplicitAtDepth(
          TIME_STEPS, MAX_TIME, DEPTH_STEPS, EFFECTIVE_DIFFUSION,
          EFFECTIVE_VELOCITY, DECAY_RATE, requested_depth,
          sample_surface_concentration_vector);
  // Test solution has the correct length
  ASSERT_EQ(TIME_STEPS + 1, depth_solution.size());
  // Test solution has the same values as the `FullyImplicit()` solution.
  for (size_t time_step = 0; time_step <= TIME_STEPS; ++time_step) {
    ASSERT_DOUBLE_EQ(solution[time_step][expected_depth_step],
                     depth_solution[time_step])
        << "Requested Depth: " << requested_depth << std::endl
        << "Expected Depth Step: " << expected_depth_step << std::endl
        << "Time Step: " << time_step;
  }
}

TEST(FullyImplicitAtDepth, RoundDown) {
  // Sample parameters
  constexpr size_t TIME_STEPS = 148;
  constexpr size_t DEPTH_STEPS = 1000;
  constexpr double MAX_TIME = 74;
  constexpr double MAX_DEPTH = 200;
  constexpr double D_STAR = 9.88011475;
  constexpr double Q_STAR = 0.3;
  constexpr double THETA_STAR = 0.2255;
  constexpr double EFFECTIVE_DIFFUSION = D_STAR / THETA_STAR;
  constexpr double EFFECTIVE_VELOCITY = Q_STAR / THETA_STAR;
  constexpr double DECAY_RATE = 0.;
  constexpr double DELTA_DEPTH = MAX_DEPTH / DEPTH_STEPS;

  const auto solution = solver::FullyImplicit<TIME_STEPS, DEPTH_STEPS>(
      MAX_TIME, MAX_DEPTH, EFFECTIVE_DIFFUSION, EFFECTIVE_VELOCITY, DECAY_RATE,
      SAMPLE_SURFACE_CONCENTRATION);

  // Convert sample surface concentrations to a vector
  std::vector<double> sample_surface_concentration_vector(
      SAMPLE_SURFACE_CONCENTRATION.size());
  std::copy(SAMPLE_SURFACE_CONCENTRATION.cbegin(),
            SAMPLE_SURFACE_CONCENTRATION.cend(),
            sample_surface_concentration_vector.begin());

  // Test at middle depth
  double requested_depth = 100.25;
  size_t expected_depth_step = 100.2 / DELTA_DEPTH;
  std::vector<double> depth_solution =
      solver::FullyImplicitAtDepth(
          TIME_STEPS, MAX_TIME, DEPTH_STEPS, EFFECTIVE_DIFFUSION,
          EFFECTIVE_VELOCITY, DECAY_RATE, requested_depth,
          sample_surface_concentration_vector);
  // Test solution against `FullyImplicit()` solution
  // Test solution has the correct length
  ASSERT_EQ(TIME_STEPS + 1, depth_solution.size());
  // Test solution has the same values
  for (size_t time_step = 0; time_step <= TIME_STEPS; ++time_step) {
    ASSERT_DOUBLE_EQ(solution[time_step][expected_depth_step],
                     depth_solution[time_step])
        << "Requested Depth: " << requested_depth << std::endl
        << "Expected Depth Step: " << expected_depth_step << std::endl
        << "Time Step: " << time_step;
  }
}

TEST(FullyImplicitAtDepth, RoundUp) {
  // Sample parameters
  constexpr size_t TIME_STEPS = 148;
  constexpr size_t DEPTH_STEPS = 1000;
  constexpr double MAX_TIME = 74;
  constexpr double MAX_DEPTH = 200;
  constexpr double D_STAR = 9.88011475;
  constexpr double Q_STAR = 0.3;
  constexpr double THETA_STAR = 0.2255;
  constexpr double EFFECTIVE_DIFFUSION = D_STAR / THETA_STAR;
  constexpr double EFFECTIVE_VELOCITY = Q_STAR / THETA_STAR;
  constexpr double DECAY_RATE = 0.;
  constexpr double DELTA_DEPTH = MAX_DEPTH / DEPTH_STEPS;

  const auto solution = solver::FullyImplicit<TIME_STEPS, DEPTH_STEPS>(
      MAX_TIME, MAX_DEPTH, EFFECTIVE_DIFFUSION, EFFECTIVE_VELOCITY, DECAY_RATE,
      SAMPLE_SURFACE_CONCENTRATION);

  // Convert sample surface concentrations to a vector
  std::vector<double> sample_surface_concentration_vector(
      SAMPLE_SURFACE_CONCENTRATION.size());
  std::copy(SAMPLE_SURFACE_CONCENTRATION.cbegin(),
            SAMPLE_SURFACE_CONCENTRATION.cend(),
            sample_surface_concentration_vector.begin());

  // Test at middle depth
  double requested_depth = 100.15;
  size_t expected_depth_step = 100.2 / DELTA_DEPTH;
  std::vector<double> depth_solution =
      solver::FullyImplicitAtDepth(
          TIME_STEPS, MAX_TIME, DEPTH_STEPS, EFFECTIVE_DIFFUSION,
          EFFECTIVE_VELOCITY, DECAY_RATE, requested_depth,
          sample_surface_concentration_vector);
  // Test solution against `FullyImplicit()` solution
  // Test solution has the correct length
  ASSERT_EQ(TIME_STEPS + 1, depth_solution.size());
  // Test solution has the same values
  for (size_t time_step = 0; time_step <= TIME_STEPS; ++time_step) {
    ASSERT_DOUBLE_EQ(solution[time_step][expected_depth_step],
                     depth_solution[time_step])
        << "Requested Depth: " << requested_depth << std::endl
        << "Expected Depth Step: " << expected_depth_step << std::endl
        << "Time Step: " << time_step;
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

  const size_t MAX_RUNTIME = 60;  // Timeout, in seconds
  if (not::testing::GTEST_FLAG(break_on_failure)) {
    signal(SIGALRM, timeout_handler);  // What to call when timer expires
    alarm(MAX_RUNTIME);                // set timer at MAX_RUNTIME seconds
  }
  return RUN_ALL_TESTS();
}
