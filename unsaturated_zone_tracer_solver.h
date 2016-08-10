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



}

namespace unsaturated_zone_tracer_solver {

namespace internal = unsaturated_zone_tracer_solver_internal;

template <size_t time_steps, size_t depth_steps>
using solution_grid = std::array<std::array<double, depth_steps>, time_steps>;

template <size_t time_steps, size_t depth_steps>
solution_grid<time_steps, depth_steps>& wieghted_time_difference(
    const unsigned int max_depth, const unsigned int max_time,
    const double effective_diffusion, const double effective_velocity,
    const std::array<double, time_steps + 1>);

}

#endif  // UNSATURATED_ZONE_TRACER_SOLVER_H
