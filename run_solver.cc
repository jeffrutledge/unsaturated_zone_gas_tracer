/**
 * \file run_solver.cc
 *
 * \author Jeffrey Rutledge jeff_rutledge@icloud.com jrutledge@usgs.gov
 *
 * \brief Uses the solver to generate a solution grid csv file.
 */

#include <iostream>
#include <fstream>
#include <limits>
#include <sstream>

#include "unsaturated_zone_tracer_solver.h"

namespace solver = unsaturated_zone_tracer_solver;
using std::cout;
using std::endl;
using std::string;

int main() {
  constexpr size_t TIME_STEPS = 148;
  constexpr size_t DEPTH_STEPS = 1000;
  constexpr double D_STAR = 9.88011475;
  constexpr double Q_STAR = 0.3;
  constexpr double THETA_STAR = 0.2255;
  constexpr double EFFECTIVE_DIFFUSION = D_STAR / THETA_STAR;
  constexpr double EFFECTIVE_VELOCITY = Q_STAR / THETA_STAR;
  constexpr double DECAY_RATE = 0.;

  // Read in surface tracer concentrations
  std::array<double, TIME_STEPS + 1> surface_tracer_concentrations;
  std::ifstream surface_tracer_csv("./cfc-11_atmospheric_concentrations.csv");
  string row;
  string entry;
  size_t entry_input_count = 0;
  getline(surface_tracer_csv, row);    // Remove column names
  while (getline(surface_tracer_csv, row) &&
         entry_input_count < TIME_STEPS + 1) {
    std::istringstream row_stream(row);
    getline(row_stream, entry, ',');  // Remove date
    getline(row_stream, entry, ',');  // Get concentration
    surface_tracer_concentrations[entry_input_count] = std::stod(entry);
    ++entry_input_count;
  }

  const auto solution = solver::CrankNicolson<TIME_STEPS, DEPTH_STEPS>(
      74, 200, EFFECTIVE_DIFFUSION, EFFECTIVE_VELOCITY, DECAY_RATE,
      surface_tracer_concentrations);
  std::ofstream solution_grid_csv;
  solution_grid_csv.open("./c++_solution_grid.csv");
  solution_grid_csv.precision(std::numeric_limits<double>::max_digits10);
  for (const auto time_step : solution) {
    auto point_iter = time_step.cbegin();
    solution_grid_csv << *point_iter;
    ++point_iter;
    for (; point_iter != time_step.cend(); ++point_iter) {
      solution_grid_csv << "," << *point_iter;
    }
    solution_grid_csv << endl;
  }
  solution_grid_csv.close();
  
  return 0;
}
