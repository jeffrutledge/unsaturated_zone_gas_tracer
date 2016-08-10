# ----- Make Variables -----

CXXFLAGS = -g -std=c++11 -Wall -Wextra -pedantic
CXX = clang++

# ----- Assignment-specific targets ---

TARGETS = run_solver

all: $(TARGETS)

# ----- Assignment-specific rules -----

unsaturated_zone_tracer_solver.o: unsaturated_zone_tracer_solver.cc
	$(CXX) $(CPPFLAGS) -c $(CXXFLAGS) unsaturated_zone_tracer_solver.cc

run_solver: run_solver.cc unsaturated_zone_tracer_solver.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

clean:
	rm -rf $(TARGETS) *.o *.dSYM documentation

doxygen:
	doxygen doxygen.config
