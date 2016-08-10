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

#include <array>

namespace unsaturated_zone_tracer_solver_internal{

/**
 * \brief Uses the Thomas Algorithm to solves a tridiagonal matrix system
 * where the matrix's lower, middle, and upper diagonal's each have a single
 * repeated value.
 *
 * Solves the system,
 * \f[ A \cdot \vec{x} = \vec{b} \f]
 *
 * \param rhs_vector The vector the matrix vector product is equal to,
 * \f$\vec{b}\f$ in the equation.
 *
 * \returns The vector, \f$\vec{x}\f$, that solves the system.
 *
 * \remark The matrix, \f$A\f$, must be diagonally dominant.
 */
template <size_t matrix_size>
std::array<double, matrix_size> thomasAlgorithimSingleValue(
    const double lower_diagonal, const double middle_diagonal,
    const double upper_diagonal,
    const std::array<double, matrix_size> rhs_vector);

}

namespace unsaturated_zone_tracer_solver {

namespace internal = unsaturated_zone_tracer_solver_internal;

/**
 * \brief An alias of a two-dimensional array used to store the solution grid.
 */
template <size_t time_steps, size_t depth_steps>
using solution_grid = std::array<std::array<double, depth_steps>, time_steps>;

/**
 * \brief Solves the PDE with the given parameters
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
 * surface boundary condition (\f$z = 0\f$) and the initial boundary condition (\f$t
 * = 0\f$).
 */
template <size_t time_steps, size_t depth_steps>
solution_grid<time_steps + 1, depth_steps + 1>& WieghtedTimeDifference(
    const unsigned int max_depth, const unsigned int max_time,
    const double effective_diffusion, const double effective_velocity,
    const std::array<double, time_steps + 1>& surface_tracer_concentrations);

}

#endif  // UNSATURATED_ZONE_TRACER_SOLVER_H
