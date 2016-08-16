/**
 * \file unsaturated_zone_tracer_solver.h
 *
 * \author Jeffrey Rutledge jeff_rutledge@icloud.com jrutledge@usgs.gov
 *
 * \brief A solver for the PDE proposed by Cook and Solomon to model the
 * transport of gas tracers through the unsaturated zone of soil.
 */

#ifndef UNSATURATED_ZONE_TRACER_SOLVER_H
#define UNSATURATED_ZONE_TRACER_SOLVER_H

#include <cassert>
#include <cmath>

#include <array>
#include <iostream>
#include <iterator>
#include <vector>

namespace unsaturated_zone_tracer_solver_internal{

/**
 * \brief Uses the Thomas Algorithm to solve a tridiagonal matrix system
 * where the matrix's lower, middle, and upper diagonal's each have a single
 * repeated value.
 *
 * Solves the system,
 * \f[ A \cdot \vec{x} = \vec{b} \f]
 *
 * \param matrix_size The dimension of the square matrix \f$A\f$. This must be
 * greater than 0.
 * \param rhs_vector The vector the matrix vector product is equal to,
 * \f$\vec{b}\f$ in the equation.
 *
 * \returns The vector, \f$\vec{x}\f$, that solves the system.
 *
 * \remark The matrix, \f$A\f$, must be diagonally dominant.
 */
template <size_t matrix_size>
std::array<double, matrix_size> ThomasAlgorithimSingleValue(
    const double lower_diagonal, const double middle_diagonal,
    const double upper_diagonal,
    const std::array<double, matrix_size> rhs_vector) {
  assert(matrix_size > 0);

  std::array<double, matrix_size> rhs_vector_prime = rhs_vector;
  std::array<double, matrix_size> middle_diagonal_prime;
  middle_diagonal_prime[0] = middle_diagonal;
  double m;
  for (size_t i = 1; i < matrix_size; ++i) {
    m = lower_diagonal / middle_diagonal_prime[i - 1];
    middle_diagonal_prime[i] = middle_diagonal - m * upper_diagonal;
    rhs_vector_prime[i] = rhs_vector_prime[i] - m * rhs_vector_prime[i - 1];
  }

  std::array<double, matrix_size> solution_vector;
  solution_vector.back() = (rhs_vector_prime.back() /
                            middle_diagonal_prime.back());
  for (size_t i = matrix_size - 1; i > 0; i--) {
    solution_vector[i - 1] = (rhs_vector_prime[i - 1] - upper_diagonal *
                              solution_vector[i]) /
                             middle_diagonal_prime[i - 1];
  }

  return solution_vector;
}

/**
 * \brief Uses the Thomas Algorithm to solves a tridiagonal matrix system
 * where the matrix's lower, middle, and upper diagonal's each have a single
 * repeated value.
 *
 * \detail Unlike the other ThomasAlgorithimSingleValue() this one is not
 * templated. Therefore, it can be used on a matrix size determined at runtime.
 *
 * Solves the system,
 * \f[ A \cdot \vec{x} = \vec{b} \f]
 *
 * \note The matrix size is determined from the size of the `rhs_vector`.
 *
 * \param matrix_size The dimension of the square matrix \f$A\f$. This must be
 * greater than 0.
 * \param rhs_vector The vector the matrix vector product is equal to,
 * \f$\vec{b}\f$ in the equation.
 *
 * \returns The vector, \f$\vec{x}\f$, that solves the system.
 *
 * \remark The matrix, \f$A\f$, must be diagonally dominant.
 */
std::vector<double> ThomasAlgorithimSingleValue(
    const double lower_diagonal, const double middle_diagonal,
    const double upper_diagonal, const std::vector<double> rhs_vector);
}

namespace unsaturated_zone_tracer_solver {

namespace internal = unsaturated_zone_tracer_solver_internal;

/// \brief An alias of a two-dimensional array used to store the solution grid.
template <size_t time_steps, size_t depth_steps>
using solution_grid = std::array<std::array<double, depth_steps>, time_steps>;

/**
 * \brief Solves the PDE using the Crank Nicolson method and returns a vector
 * of the concentrations at the depth step closest to the requested depth.
 *
 * \detail Unlike CrankNicolson() this function is not templated on `time_steps`
 * and `depth_steps`. This means these may be specified at run time instead of
 * at compile time. Another difference is the max depth is always 200m.
 *
 * [This pdf](./solver_method.pdf) details the numerical method used, as well as
 * discusses its accuracy and stability.
 *
 * \pre `surface_tracer_concentrations` must contain atleast `time_steps` + 1
 * elements. One extra for the initial boundary and then one for each time step
 * the solver will take.
 *
 * \param time_steps Number of time steps to take. This can be calculated as the
 * product of `max_time` and the time step size.
 * \param max_time The final time to step to. Starting at \f$t = 0\f$, meaning
 * if the data starts at 1940 and ends at 2000.5 the max_time would be \f$2014 -
 * 1940 = 74\f$.
 * \param depth_steps Number of depth steps to take.
 * \param effective_diffusion The constant coefficient of the second order
 * partial of concentration with respect to position.
 * \param effective_velocity The constant coefficient of the first order
 * partial of concentration with respect to position.
 * \param decay_constant The rate at which the species decays.
 * \param surface_tracer_concentration The concentration of the tracer gas at
 * the surface. This should contain one more element than the number of time
 * steps for a boundary condition and it should have one value for each time
 * step.
 *
 * \return A vector containing the concentrations at every time step, including
 * the initial boundary, at the depth step closest to the requested depth. It
 * will be of length `time_steps` + 1.
 */
std::vector<double> CrankNicolsonAtDepth(
    const size_t time_steps, const double max_time, const size_t depth_steps,
    const double effective_diffusion, const double effective_velocity,
    const double decay_rate, const double requested_depth,
    const std::vector<double>& surface_tracer_concentrations);

/**
 * \brief Solves the PDE using a Crank Nicolson finite difference scheme.
 *
 * \detail Time and depth will step from 0 to the max parameter in equal sized
 * increments. The depth step after the `max_depth` is bounded at a
 * concentration of 0.
 *
 * [This pdf](./solver_method.pdf) details the method used, as well as
 * discusses its accuracy and stability.
 *
 * \param time_steps Number of time steps to take. This can be calculated as the
 * product of `max_time` and the time step size.
 * \param depth_steps Number of depth steps to take.
 * \param max_depth The final depth to step to.
 * \param max_time The final time to step to. Starting at \f$t = 0\f$, meaning
 * if the data starts at 1940 and ends at 2000.5 the max_time would be \f$2014 -
 * 1940 = 74\f$.
 * \param effective_diffusion The constant coefficient of the second order
 * partial of concentration with respect to position.
 * \param effective_velocity The constant coefficient of the first order
 * partial of concentration with respect to position.
 * \param decay_constant The rate at which the species decays.
 * \param surface_tracer_concentration The concentration of the tracer gas at
 * the surface. This must contain one more element than the number of time
 * steps for a boundary condition and then it should have one value for each
 * time step.
 *
 * \return A two dimensional array of the solution grid containing all the
 * surface boundary condition (\f$z = 0\f$) and the initial boundary condition
 * (\f$t = 0\f$).
 */
template <size_t time_steps, size_t depth_steps>
solution_grid<time_steps + 1, depth_steps + 1> CrankNicolson(
    const double max_time, const double max_depth,
    const double effective_diffusion, const double effective_velocity,
    const double decay_rate,
    const std::array<double, time_steps + 1>& surface_tracer_concentrations) {
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
  // Initialize the solution grid with boundary condition at t = 0
  solution_grid<time_steps + 1, depth_steps + 1> solution;
  solution[0].fill(0.);
  solution[0][0] = surface_tracer_concentrations[0];

  std::array<double, depth_steps> previous_time_vector;
  for (size_t time_step = 1; time_step < time_steps + 1; ++time_step) {
    // Add surface concentration boundary to solution
    solution[time_step][0] = surface_tracer_concentrations[time_step];
    // Calculate the RHS vector from the previous time step. This calculates the
    // product of the tridiagonal matrix of the previous time step and the
    // previous time step solution vector, and adds the boundary offset.
    previous_time_vector[0] =
        (previous_time_middle_diagonal * solution[time_step - 1][1] +
         previous_time_upper_diagonal * solution[time_step - 1][2] +
         // Now add boundary condition offsets
         previous_time_lower_diagonal * solution[time_step - 1][0] -
         current_time_lower_diagonal * solution[time_step][0]);
    for (size_t row = 1; row < previous_time_vector.size() - 1; ++row) {
      previous_time_vector[row] =
          (previous_time_lower_diagonal * solution[time_step - 1][row] +
           previous_time_middle_diagonal * solution[time_step - 1][row + 1] +
           previous_time_upper_diagonal * solution[time_step - 1][row + 2]);
    }
    // Boundary offset here is 0 because boundary at max_depth is 0.
    previous_time_vector.back() = 
        (previous_time_lower_diagonal *
         solution[time_step - 1][depth_steps - 1] +
         previous_time_middle_diagonal * solution[time_step - 1][depth_steps]);

    // Calculate solution for this time step and insert it into solution
    const auto insert_start = std::next(solution[time_step].begin(), 1);
    std::array<double, depth_steps> current_time_solution =
        internal::ThomasAlgorithimSingleValue<depth_steps>(
            current_time_lower_diagonal, current_time_middle_diagonal,
            current_time_upper_diagonal, previous_time_vector);
    std::copy(current_time_solution.begin(), current_time_solution.end(),
              insert_start);
  }
  return solution;
}

/**
 * \brief Solves the PDE using a fully implicit finite difference scheme and
 * returns a vector of the concentrations at the depth step closest to the
 * requested depth.
 *
 * \detail Unlike FullyImplicit() this function is not templated on `time_steps`
 * and `depth_steps`. This means these may be specified at run time instead of
 * at compile time. Another difference is the max depth is always 200m.
 *
 * [This pdf](./solver_method.pdf) details the numerical method used, as well as
 * discusses its accuracy and stability.
 *
 * \pre `surface_tracer_concentrations` must contain atleast `time_steps` + 1
 * elements. One extra for the initial boundary and then one for each time step
 * the solver will take.
 *
 * \param time_steps Number of time steps to take. This can be calculated as the
 * product of `max_time` and the time step size.
 * \param max_time The final time to step to. Starting at \f$t = 0\f$, meaning
 * if the data starts at 1940 and ends at 2000.5 the max_time would be \f$2014 -
 * 1940 = 74\f$.
 * \param depth_steps Number of depth steps to take.
 * \param effective_diffusion The constant coefficient of the second order
 * partial of concentration with respect to position.
 * \param effective_velocity The constant coefficient of the first order
 * partial of concentration with respect to position.
 * \param decay_constant The rate at which the species decays.
 * \param surface_tracer_concentration The concentration of the tracer gas at
 * the surface. This should contain one more element than the number of time
 * steps for a boundary condition and it should have one value for each time
 * step.
 *
 * \return A vector containing the concentrations at every time step, including
 * the initial boundary, at the depth step closest to the requested depth. It
 * will be of length `time_steps` + 1.
 */
std::vector<double> FullyImplicitAtDepth(
    const size_t time_steps, const double max_time, const size_t depth_steps,
    const double effective_diffusion, const double effective_velocity,
    const double decay_rate, const double requested_depth,
    const std::vector<double>& surface_tracer_concentrations);

/**
 * \brief Solves the PDE with the given parameters using a fully implicit finite
 * difference scheme.
 *
 * The standard forward difference time approximation,
 * \f[
 *  \frac{\partial u}{\partial t} \approx
 *    \frac{u_{m+1,\ j} - u_{m,\ j}}{\Delta t}.
 * \f]
 *
 * \detail Time and depth will step from 0 to the max parameter in equal sized
 * increments. The depth step after the `max_depth` is bounded at a
 * concentration of 0.
 *
 * [This pdf](./solver_method.pdf) details the method used, as well as
 * discusses its accuracy and stability.
 *
 * \param time_steps Number of time steps to take. This can be calculated as the
 * product of `max_time` and the time step size.
 * \param depth_steps Number of depth steps to take.
 * \param max_depth The final depth to step to.
 * \param max_time The final time to step to. Starting at \f$t = 0\f$, meaning
 * if the data starts at 1940 and ends at 2000.5 the max_time would be \f$2014 -
 * 1940 = 74\f$.
 * \param effective_diffusion The constant coefficient of the second order
 * partial of concentration with respect to position.
 * \param effective_velocity The constant coefficient of the first order
 * partial of concentration with respect to position.
 * \param decay_constant The rate at which the species decays.
 * \param surface_tracer_concentration The concentration of the tracer gas at
 * the surface. This must contain one more element than the number of time
 * steps for a boundary condition and then it should have one value for each
 * time step.
 *
 * \return A two dimensional array of the solution grid containing all the
 * surface boundary condition (\f$z = 0\f$) and the initial boundary condition
 * (\f$t = 0\f$).
 */
template <size_t time_steps, size_t depth_steps>
solution_grid<time_steps + 1, depth_steps + 1> FullyImplicit(
    const double max_time, const double max_depth,
    const double effective_diffusion, const double effective_velocity,
    const double decay_rate,
    const std::array<double, time_steps + 1>& surface_tracer_concentrations) {
  const double delta_time = max_time / time_steps;
  const double delta_depth = max_depth / depth_steps;
  // Construct factors in the iterative equation obtained from the finite
  // difference method. These parameters are described in the
  // solver_method.tex file.
  // TODO Update the solver_methods.tex with these parameters.
  const double alpha = (effective_diffusion * delta_time) /
      (std::pow(delta_depth, 2));
  const double beta = (effective_velocity * delta_time) / (2 * delta_depth);
  // Construct the diagonal entries of the tridiagonal matrices.
  // For the current time step:
  const double current_time_lower_diagonal = -alpha - beta;
  const double current_time_middle_diagonal = 1 + decay_rate + 2 * alpha;
  const double current_time_upper_diagonal = -alpha + beta;
  // Initialize the solution grid with boundary condition at t = 0
  solution_grid<time_steps + 1, depth_steps + 1> solution;
  solution[0].fill(0.);
  solution[0][0] = surface_tracer_concentrations[0];

  std::array<double, depth_steps> previous_time_vector;
  for (size_t time_step = 1; time_step < time_steps + 1; ++time_step) {
    // Add surface concentration boundary to solution.
    solution[time_step][0] = surface_tracer_concentrations[time_step];
    // Calculate the RHS vector from the previous time step. This is simply the
    // previous time step solution, plus some boundary offsets.
    std::copy(std::next(solution[time_step - 1].begin(), 1),
              std::prev(solution[time_step - 1].end(), 1),
              previous_time_vector.begin());
     // Now add boundary condition offsets.
     // Note, the boundary offset at max depth is 0 because boundary at
     // max_depth is 0.
    previous_time_vector[0] -= (current_time_lower_diagonal *
                                solution[time_step][0]);

    // Calculate the solution for this time step and insert it into the solution
    // grid.
    const auto insert_start = std::next(solution[time_step].begin(), 1);
    std::array<double, depth_steps> current_time_solution =
        internal::ThomasAlgorithimSingleValue<depth_steps>(
            current_time_lower_diagonal, current_time_middle_diagonal,
            current_time_upper_diagonal, previous_time_vector);
    std::copy(current_time_solution.begin(), current_time_solution.end(),
              insert_start);
  }
  return solution;
}

/**
 * \brief Solves the PDE with the given parameters using the weighted forward
 * time approximation.
 *
 * The weighted forward time difference approximation used is,
 * \f[
 *  \frac{\partial u}{\partial t} \approx
 *    \frac{u_{m+1,\ j+1} - u_{m,\ j+1}}{6\Delta t} +
 *    2\frac{u_{m+1,\ j} - u_{m,\ j}}{3\Delta t} +
 *    \frac{u_{m+1,\ j-1} - u_{m,\ j-1}}{6\Delta t}.
 * \f]
 *
 * \detail Time and depth will step from 0 to the max parameter in equal sized
 * increments. The depth step after the `max_depth` is bounded at a
 * concentration of 0.
 *
 * \param time_steps Number of time steps to take.
 * \param depth_steps Number of depth steps to take.
 * \param max_depth The final depth to step to.
 * \param max_time The final time to step to.
 * \param effective_diffusion The constant coefficient of the second order
 * partial of concentration with respect to position.
 * \param effective_velocity The constant coefficient of the first order
 * partial of concentration with respect to position.
 * \param surface_tracer_concentration The concentration of the tracer gas at
 * the surface. This should contain one more element than the number of time
 * steps for a boundary condition and it should have one value for each time
 * step.
 *
 * \return A two dimensional array of the solution grid containing all the
 * surface boundary condition (\f$z = 0\f$) and the initial boundary condition
 * (\f$t = 0\f$).
 */
template <size_t time_steps, size_t depth_steps>
solution_grid<time_steps + 1, depth_steps + 1> WieghtedTimeDifference(
    const double max_time, const double max_depth,
    const double effective_diffusion, const double effective_velocity,
    const std::array<double, time_steps + 1>& surface_tracer_concentrations) {
  const double delta_time = max_time / time_steps;
  const double delta_depth = max_depth / depth_steps;
  // Construct factors in the iterative equation obtained from the finite
  // difference method.
  const double time_approx_factor = std::pow(delta_depth, 2) / delta_time;
  // Note: These remaining factors are for the current time step side of the
  // equation, so they will be negated when used in calculations on the previous
  // time step side.
  const double previous_depth_factor = (-6 * effective_diffusion - 3 *
                                        effective_velocity * delta_depth);
  const double current_depth_factor = 12 * effective_diffusion;
  const double next_depth_factor = (-6 * effective_diffusion + 3 *
                                    effective_velocity * delta_depth);
  // Construct the diagonal entries of the tridiagonal matrices.
  // For the current time step:
  const double current_time_lower_diagonal = (2 * time_approx_factor +
                                              previous_depth_factor);
  const double current_time_middle_diagonal = (8 * time_approx_factor +
                                               current_depth_factor);
  const double current_time_upper_diagonal = (2 * time_approx_factor +
                                              next_depth_factor);
  // For the previous time step:
  const double previous_time_lower_diagonal = (2 * time_approx_factor -
                                               previous_depth_factor);
  const double previous_time_middle_diagonal = (8 * time_approx_factor -
                                                current_depth_factor);
  const double previous_time_upper_diagonal = (2 * time_approx_factor -
                                               next_depth_factor);
  // Initialize the solution grid with boundary condition at t = 0
  solution_grid<time_steps + 1, depth_steps + 1> solution;
  solution[0].fill(0.);
  solution[0][0] = surface_tracer_concentrations[0];

  std::array<double, depth_steps> previous_time_vector;
  for (size_t time_step = 1; time_step < time_steps + 1; ++time_step) {
    // Add surface concentration boundary to solution
    solution[time_step][0] = surface_tracer_concentrations[time_step];
    // Calculate the RHS vector from the previous time step. This calculates the
    // product of the tridiagonal matrix of the previous time step and the
    // previous time step solution vector, and adds the boundary offset.
    previous_time_vector[0] =
        (previous_time_middle_diagonal * solution[time_step - 1][1] +
         previous_time_upper_diagonal * solution[time_step - 1][2] +
         // Now add boundary condition offsets
         previous_time_lower_diagonal * solution[time_step - 1][0] -
         current_time_lower_diagonal * solution[time_step][0]);
    for (size_t row = 1; row < previous_time_vector.size() - 1; ++row) {
      previous_time_vector[row] =
          (previous_time_lower_diagonal * solution[time_step - 1][row] +
           previous_time_middle_diagonal * solution[time_step - 1][row + 1] +
           previous_time_upper_diagonal * solution[time_step - 1][row + 2]);
    }
    // Boundary offset here is 0 because boundary at max_depth is 0.
    previous_time_vector.back() = 
        (previous_time_lower_diagonal *
         solution[time_step - 1][depth_steps - 1] +
         previous_time_middle_diagonal * solution[time_step - 1][depth_steps]);

    // Calculate solution for this time step and insert it into solution
    const auto insert_start = std::next(solution[time_step].begin(), 1);
    std::array<double, depth_steps> current_time_solution =
        internal::ThomasAlgorithimSingleValue<depth_steps>(
            current_time_lower_diagonal, current_time_middle_diagonal,
            current_time_upper_diagonal, previous_time_vector);
    std::copy(current_time_solution.begin(), current_time_solution.end(),
              insert_start);
  }
  return solution;
}

}

#endif  // UNSATURATED_ZONE_TRACER_SOLVER_H
