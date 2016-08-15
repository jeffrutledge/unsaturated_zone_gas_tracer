# ----- Make Variables -----

CXXFLAGS = -g -std=c++11 -Wall -Wextra -pedantic
CXX = clang++

# ----- Google Test Setup -----

GTEST_DIR = /usr/local/include/googletest
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/gtest*.h \
		$(GTEST_DIR)/include/gtest/internal/*.h
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)
CPPFLAGS = -isystem $(GTEST_DIR)/include

# ----- Assignment-specific targets ---

TARGETS = run_solver test_unsaturated_zone_tracer_solver

all: $(TARGETS)

# ----- Google Test Rules -----

gtest-all.o: $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
		-pthread $(GTEST_DIR)/src/gtest-all.cc

# ----- Assignment-specific rules -----

unsaturated_zone_tracer_solver.o: unsaturated_zone_tracer_solver.cc
	$(CXX) $(CPPFLAGS) -c $(CXXFLAGS) unsaturated_zone_tracer_solver.cc

run_solver: run_solver.cc unsaturated_zone_tracer_solver.o unsaturated_zone_tracer_solver.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

test_unsaturated_zone_tracer_solver: test_unsaturated_zone_tracer_solver.cc unsaturated_zone_tracer_solver.o gtest-all.o unsaturated_zone_tracer_solver.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< unsaturated_zone_tracer_solver.o -o $@ -pthread gtest-all.o

test: test_unsaturated_zone_tracer_solver
	./test_unsaturated_zone_tracer_solver

clean:
	rm -rf $(TARGETS) *.o *.dSYM documentation gtest.a gtest_main.a

doxygen:
	doxygen doxygen.config
