/**
 * \file unsaturated_zone_tracer_solver.cc
 *
 * \author Jeffrey Rutledge jeff_rutledge@icloud.com jrutledge@usgs.gov
 *
 * \brief A solver for the PDE proposed by Cook and Solomon to model the
 * transport of gas tracers through the unsaturated zone of soil.
 */
#include "unsaturated_zone_tracer_solver.h"

#include <math.h>

namespace unsaturated_zone_tracer_solver_internal {

std::vector<double> ThomasAlgorithimSingleValue(
    const double lower_diagonal, const double middle_diagonal,
    const double upper_diagonal,
    const std::vector<double> rhs_vector) {
  const size_t matrix_size = rhs_vector.size();
  assert(matrix_size > 0);

  std::vector<double> rhs_vector_prime = rhs_vector;
  std::vector<double> middle_diagonal_prime(matrix_size);
  middle_diagonal_prime[0] = middle_diagonal;
  double m;
  for (size_t i = 1; i < matrix_size; ++i) {
    m = lower_diagonal / middle_diagonal_prime[i - 1];
    middle_diagonal_prime[i] = middle_diagonal - m * upper_diagonal;
    rhs_vector_prime[i] = rhs_vector_prime[i] - m * rhs_vector_prime[i - 1];
  }

  std::vector<double> solution_vector(matrix_size);
  solution_vector.back() = (rhs_vector_prime.back() /
                            middle_diagonal_prime.back());
  for (size_t i = matrix_size - 1; i > 0; i--) {
    solution_vector[i - 1] = (rhs_vector_prime[i - 1] - upper_diagonal *
                              solution_vector[i]) /
                             middle_diagonal_prime[i - 1];
  }

  return solution_vector;
}

}

namespace unsaturated_zone_tracer_solver {

std::vector<double> CrankNicolsonAtDepth(
    const size_t time_steps, const double max_time, const size_t depth_steps,
    const double effective_diffusion, const double effective_velocity,
    const double decay_rate, const double requested_depth,
    const std::vector<double>& surface_tracer_concentrations) {
  const double max_depth = 200.;
  const double delta_time = max_time / time_steps;
  const double delta_depth = max_depth / depth_steps;
  // Construct factors in the iterative equation obtained from the finite
  // difference method.
  // Note: These remaining factors are for the current time step side of the
  // equation, so they will be negated when used in calculations on the previous
  // time step side.
  const double alpha = (effective_diffusion * delta_time) /
      (2 * std::pow(delta_depth, 2));
  const double beta = (effective_velocity * delta_time) / (4 * delta_depth);
  // Construct the diagonal entries of the tridiagonal matrices.
  // For the current time step:
  const double current_time_lower_diagonal = -alpha - beta;
  const double current_time_middle_diagonal = 1 + decay_rate + 2 * alpha;
  const double current_time_upper_diagonal = -alpha + beta;
  // For the previous time step:
  const double previous_time_lower_diagonal = alpha + beta;
  const double previous_time_middle_diagonal = 1 - 2 * alpha;
  const double previous_time_upper_diagonal = alpha - beta;

  // Initialize the previous time solution with boundary condition at t = 0
  std::vector<double> previous_time_solution(depth_steps + 1, 0.);
  previous_time_solution[0] = surface_tracer_concentrations[0];
  // Find the depth step closest to the requested depth
  size_t requested_depth_step = round(requested_depth / delta_depth);
  std::vector<double> solution_at_requested_depth(time_steps + 1);
  // Add initial value to solution_at_requested_depth;
  solution_at_requested_depth[0] = previous_time_solution[requested_depth_step];

  std::vector<double> previous_time_vector(depth_steps);
  for (size_t time_step = 1; time_step < time_steps + 1; ++time_step) {
    // Calculate the RHS vector from the previous time step. This calculates the
    // product of the tridiagonal matrix of the previous time step and the
    // previous time step solution vector, and adds the boundary offset.
    previous_time_vector[0] =
        (previous_time_middle_diagonal * previous_time_solution[1] +
         previous_time_upper_diagonal * previous_time_solution[2] +
         // Now add boundary condition offsets
         previous_time_lower_diagonal * previous_time_solution[0] -
         current_time_lower_diagonal *
         surface_tracer_concentrations[time_step]);
    for (size_t row = 1; row < previous_time_vector.size() - 1; ++row) {
      previous_time_vector[row] =
          (previous_time_lower_diagonal * previous_time_solution[row] +
           previous_time_middle_diagonal * previous_time_solution[row + 1] +
           previous_time_upper_diagonal * previous_time_solution[row + 2]);
    }
    // Boundary offset here is 0 because boundary at max_depth is 0.
    previous_time_vector.back() = 
        (previous_time_lower_diagonal *
         previous_time_solution[depth_steps - 1] +
         previous_time_middle_diagonal * previous_time_solution[depth_steps]);

    // Calculate solution for this time step
    std::vector<double> current_time_solution =
        internal::ThomasAlgorithimSingleValue(
            current_time_lower_diagonal, current_time_middle_diagonal,
            current_time_upper_diagonal, previous_time_vector);
    // Insert the solution into the previous_time_solution
    const auto insert_start = std::next(previous_time_solution.begin(), 1);
    std::copy(current_time_solution.begin(), current_time_solution.end(),
              insert_start);
    // Set surface boundary value of previous_time_solution
    previous_time_solution[0] = surface_tracer_concentrations[time_step];
    // Add the relevant value to the wanted_depth_solution
    solution_at_requested_depth[time_step] =
        previous_time_solution[requested_depth_step];
  }
  return solution_at_requested_depth;
}

std::vector<double> FullyImplicitAtDepth(
    const size_t time_steps, const double max_time, const size_t depth_steps,
    const double effective_diffusion, const double effective_velocity,
    const double decay_rate, const double requested_depth,
    const std::vector<double>& surface_tracer_concentrations) {
  const double max_depth = 200.;
  const double delta_time = max_time / time_steps;
  const double delta_depth = max_depth / depth_steps;
  // Construct factors in the iterative equation obtained from the finite
  // difference method. These parameters are described in the
  // solver_method.tex file.
  const double alpha = (effective_diffusion * delta_time) /
      (std::pow(delta_depth, 2));
  const double beta = (effective_velocity * delta_time) / (2 * delta_depth);
  // Construct the diagonal entries of the tridiagonal matrix.
  const double current_time_lower_diagonal = -alpha - beta;
  const double current_time_middle_diagonal = 1 + decay_rate + 2 * alpha;
  const double current_time_upper_diagonal = -alpha + beta;

  // Initialize the previous time solution with boundary condition at t = 0
  std::vector<double> previous_time_solution(depth_steps + 1, 0.);
  previous_time_solution[0] = surface_tracer_concentrations[0];
  // Find the depth step closest to the requested depth
  size_t requested_depth_step = round(requested_depth / delta_depth);
  std::vector<double> solution_at_requested_depth(time_steps + 1);
  // Add initial value to solution_at_requested_depth;
  solution_at_requested_depth[0] = previous_time_solution[requested_depth_step];
  std::vector<double> previous_time_vector(depth_steps);
  for (size_t time_step = 1; time_step < time_steps + 1; ++time_step) {
    // Calculate the RHS vector from the previous time step. This is simply the
    // previous time step solution, plus some boundary offsets.
    std::copy(std::next(previous_time_solution.begin(), 1),
              std::prev(previous_time_solution.end(), 1),
              previous_time_vector.begin());
     // Now add boundary condition offsets.
     // Note, the boundary offset at max depth is 0 because boundary at
     // max_depth is 0.
    previous_time_vector[0] -= (current_time_lower_diagonal *
                                surface_tracer_concentrations[time_step]);

    // Calculate solution for this time step
    std::vector<double> current_time_solution =
        internal::ThomasAlgorithimSingleValue(
            current_time_lower_diagonal, current_time_middle_diagonal,
            current_time_upper_diagonal, previous_time_vector);
    // Insert the solution into the previous_time_solution
    const auto insert_start = std::next(previous_time_solution.begin(), 1);
    std::copy(current_time_solution.begin(), current_time_solution.end(),
              insert_start);
    // Set surface boundary value of previous_time_solution
    previous_time_solution[0] = surface_tracer_concentrations[time_step];
    // Add the relevant value to the wanted_depth_solution
    solution_at_requested_depth[time_step] =
        previous_time_solution[requested_depth_step];
  }
  return solution_at_requested_depth;
}

}
