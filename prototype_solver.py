"""
Planning for implementation of a numerical solver using the Crank
Nicolson method. For the diffusion/advection PDE of gas tracers through
the unsaturated zone.

Author: Jeffrey Rutledge (jrutledge@usgs.gov, jeff_rutledge@icloud.com)
"""

import numpy as np
from scipy.sparse import diags
from matplotlib import pyplot as plt

import timeit

class TridiagonalMatrix(object):
    def __init__(self, size, upper=None, middle=None, lower=None):
        self.size_ = size
        if upper:
            assert len(upper) == size - 1, (
                'upper does not have correct dimensions.'
                ' It is {} not {}'.format(len(upper), size - 1))
            self.upper_ = np.array(upper)
        else:
            self.upper_ = np.zeros(size - 1)

        if middle:
            assert len(middle) == size, (
                'middle does not have correct dimensions.'
                ' It is {} not {}'.format(len(middle), size))
            self.middle_ = np.array(middle)
        else:
            self.middle_ = np.zeros(size)

        if lower:
            assert len(lower) == size - 1, (
                'lower does not have correct dimensions.'
                ' It is {} not {}'.format(len(lower), size - 1))
            self.lower_ = np.array(lower)
        else:
            self.lower_ = np.zeros(size - 1)

    def get(self, row, column):
        assert row >= 0 and row < self.size_, ('Index out of bounds.'
                                              'Row {}'.format(row))
        assert column >= 0 and column < self.size_, ('Index out of bounds.'
                                                    'Column {}'.format(column))

        if row == column:
            return self.middle_[row]
        elif row == column - 1:
            return self.upper_[row]
        elif row == column + 1:
            return self.lower_[column]
        else:
            return 0.

    def set(self, row, column, value):
        assert row >= 0 and row < self.size_, ('Index out of bounds.'
                                              'Row {}'.format(row))
        assert column >= 0 and column < self.size_, ('Index out of bounds.'
                                                    'Column {}'.format(column))

        if row == column:
            self.middle_[row] = value
        elif row == column - 1:
            self.upper_[row] = value
        elif row == column + 1:
            self.lower_[column] = value
        else:
            assert False, 'Index not settable ({}, {})'.format(row, column)

    def size(self):
        return self.size_

    def __str__(self):
        output = ''
        for row in range(self.size_):
            for column in range(self.size_):
                output += '{:>8.2}'.format(self.get(row, column))
            output += '\n'
        return output


def crout_factorization(a, b):
    """
    Uses the Crout Factorization algorithm to solve the linear system
        a . x = b
    where a is a tridiagonal matrix that is positive definite, or
    strictly diagonally dominant, b is a column vector of the same
    dimension as a, for the column vector x.

    Reference
    ---------
    Richard L. Burden and J. Douglas Faires, Numerical Analysis 8th Edition,
    pp. 408-409, 2005.
    """
    n = a.size()
    # Initialize the matrices the input will be factored into: l and u
    l = TridiagonalMatrix(n)
    u = TridiagonalMatrix(n, middle=[1.] * n)

    z = np.zeros(n)

    l.set(0, 0, a.get(0, 0))
    u.set(0, 1, a.get(0, 1) / l.get(0, 0))
    z[0] = b[0] / l.get(0, 0)

    for i in range(1, n - 1):
        l.set(i, i - 1, a.get(i, i - 1))
        l.set(i, i, a.get(i, i) - l.get(i, i - 1) * u.get(i - 1, i))
        u.set(i, i + 1, a.get(i, i + 1) / l.get(i, i))
        z[i] = (b[i] - l.get(i, i - 1) * z[i - 1]) / l.get(i, i)

    l.set(n - 1, n - 2, a.get(n - 1, n - 2))
    l.set(n - 1, n - 1, a.get(n - 1, n - 1) - l.get(n - 1, n - 2) *
          u.get(n - 2, n - 1))
    z[n - 1] = ((b[n - 1] - l.get(n - 1, n - 2) * z[n - 2]) /
                l.get(n - 1, n - 1))

    output = np.zeros(n)
    output[n - 1] = z[n - 1]
    for i in range(n - 2, -1, -1):
        output[i] = z[i] - u.get(i, i + 1) * output[i + 1]
    return output


def crank_nicolson(max_depth, max_time, depth_steps, time_steps,
                   effective_diffusion_constant, effective_velocity_constant,
                   surface_tracer_concentrations):
    """
    Uses the Crank Nicolson method to numerically approximate the PDE,
        du/dt = ED (d2u/dz2) - EV (du/dz)
    where the `d`s are partials, ED is the effective diffusion constant,
    EV is the effective velocity constant, z is depth, t is time, and u
    is the gas concentration.

    The surface tracer concentration is the bound for u(t, 0) where
    u(t, z). This must have a point for each step in the simulation and one for
    corner bound at u(0, 0).

    The boundaries along u(0, z) and u(t, max_depth) are assumed to be
    u = 0.

    Method
    ------
    The PDE is approximated using the central difference approximation
    for the second and first order partials of u with respect to z. Then
    using the Crank Nicolson method by using the forward difference
    approximation for the partial of u with respect to t averaging the
    previous approximations at the current and next time step.

    These approximations give a linear system with a tridiagonal matrix
    that is solved using the Crout Factorization algorithm.
    """
    assert len(surface_tracer_concentrations) == time_steps + 1, \
        'Surface boundary has length {}, not the expected {}'.format(
            len(surface_tracer_concentrations), time_steps + 1)
    delta_time = float(max_time) / time_steps
    delta_depth = float(max_depth) / depth_steps

    # The factors in iterative equation provided by the finite difference
    # approximations
    step_factor = np.power(delta_depth, 2) / delta_time
    # Note: These factors are for the current time step side of the equation
    # so the will be negated for calculations of the previous time step
    previous_depth_factor = (-6 * effective_diffusion_constant - 3 *
                             effective_velocity_constant * delta_depth)
    current_depth_factor = 12 * effective_diffusion_constant
    next_depth_factor = (-6 * effective_diffusion_constant + 3 *
                         effective_velocity_constant * delta_depth)

    # Construct the tridiagonal matrix for the current time step side
    current_time_upper_diagonal = next_depth_factor
    current_time_middle_diagonal = 12 * step_factor + current_depth_factor
    current_time_lower_diagonal = previous_depth_factor
    current_time_matrix = TridiagonalMatrix(
        depth_steps, [current_time_upper_diagonal] * (depth_steps - 1),
        [current_time_middle_diagonal] * (depth_steps),
        [current_time_lower_diagonal] * (depth_steps - 1))

    # Construct the tridiagonal matrix for the previous time step side
    previous_time_upper_diagonal = -next_depth_factor
    previous_time_middle_diagonal = 12 * step_factor - current_depth_factor
    previous_time_lower_diagonal = -previous_depth_factor
    previous_time_matrix = \
        (np.diag([previous_time_upper_diagonal] * (depth_steps - 1), 1) +
         np.diag([previous_time_middle_diagonal] * (depth_steps), 0) +
         np.diag([previous_time_lower_diagonal] * (depth_steps - 1), -1))

    # Initialize the solution grid matrix. This matrix will contain all the
    # solution values, including the boundaries.
    solution_grid = np.empty((time_steps + 1, depth_steps + 1))
    # Fill with boundaries
    solution_grid[0,:] = 0
    solution_grid[:,0] = surface_tracer_concentrations
    for time_step in range(1, time_steps + 1):
        # Slice the boundary conditions from previous time step solution
        b = previous_time_matrix.dot(solution_grid[time_step - 1, 1:])

        # Add the boundary offset to b
        # (the offset for b[-1] is 0 because the boundary at max depth is 0)
        b[0] += (solution_grid[time_step - 1, 0] *
                 previous_time_lower_diagonal -
                 solution_grid[time_step, 0] *
                 current_time_lower_diagonal)

        # Calculate the current time step's solution
        # and insert it in between boundaries in solution matrix
        solution_grid[time_step, 1:] = crout_factorization(
            current_time_matrix, b)

    return solution_grid


def wieghted_time(max_depth, max_time, depth_steps, time_steps,
            effective_diffusion_constant, effective_velocity_constant,
            surface_tracer_concentrations):
    """
    Uses the Crank Nicolson method to numerically approximate the PDE,
        du/dt = ED (d2u/dz2) - EV (du/dz)
    where the `d`s are partials, ED is the effective diffusion constant,
    EV is the effective velocity constant, z is depth, t is time, and u
    is the gas concentration.

    The surface tracer concentration is the bound for u(t, 0) where
    u(t, z). This must have a point for each step in the simulation and
    one for corner bound at u(0, 0).

    The boundaries along u(0, z) and u(t, max_depth + one_step) are
    assumed to be u = 0.

    Method
    ------
    The PDE is approximated using the central difference approximation
    for the second and first order partials of u with respect to z. Then
    using the Crank Nicolson method by using a weighted forward

    These approximations give a linear system with a tridiagonal matrix
    that is solved using the Crout Factorization algorithm.
    """
    assert len(surface_tracer_concentrations) == time_steps + 1, \
        'Surface boundary has length {}, not the expected {}'.format(
            len(surface_tracer_concentrations), time_steps + 1)
    delta_time = float(max_time) / time_steps
    delta_depth = float(max_depth) / depth_steps

    # The factors in iterative equation provided by the finite difference
    # approximations
    step_factor = np.power(delta_depth, 2) / delta_time
    # Note: These factors are for the current time step side of the equation
    # so the will be negated for calculations of the previous time step
    previous_depth_factor = (-6 * effective_diffusion_constant - 3 *
                             effective_velocity_constant * delta_depth)
    current_depth_factor = 12 * effective_diffusion_constant
    next_depth_factor = (-6 * effective_diffusion_constant + 3 *
                         effective_velocity_constant * delta_depth)

    # Construct the tridiagonal matrix for the current time step side
    current_time_upper_diagonal = 2 * step_factor + next_depth_factor
    current_time_middle_diagonal = 8 * step_factor + current_depth_factor
    current_time_lower_diagonal = 2 * step_factor + previous_depth_factor
    current_time_matrix = TridiagonalMatrix(
        depth_steps, [current_time_upper_diagonal] * (depth_steps - 1),
        [current_time_middle_diagonal] * (depth_steps),
        [current_time_lower_diagonal] * (depth_steps - 1))

    # Construct the tridiagonal matrix for the previous time step side
    previous_time_upper_diagonal = 2 * step_factor - next_depth_factor
    previous_time_middle_diagonal = 8 * step_factor - current_depth_factor
    previous_time_lower_diagonal = 2 * step_factor - previous_depth_factor
    previous_time_matrix = \
        (np.diag([previous_time_upper_diagonal] * (depth_steps - 1), 1) +
         np.diag([previous_time_middle_diagonal] * (depth_steps), 0) +
         np.diag([previous_time_lower_diagonal] * (depth_steps - 1), -1))

    # # TESTING
    # current_time_matrix.set(
    #         current_time_matrix.size() - 1, current_time_matrix.size() - 1,
    #          4 * step_factor - previous_depth_factor)

    # previous_time_matrix[-1, -1] = 4 * step_factor + previous_depth_factor
    # ###############

    # Initialize the solution grid matrix. This matrix will contain all the
    # solution values, including the boundaries.
    solution_grid = np.empty((time_steps + 1, depth_steps + 1))
    # Fill with boundaries
    solution_grid[0,:] = 0
    solution_grid[:,0] = surface_tracer_concentrations
    for time_step in range(1, time_steps + 1):
        # Slice the boundary conditions from previous time step solution
        b = previous_time_matrix.dot(solution_grid[time_step - 1, 1:])

        # Add the boundary offset to b
        # (the offset for b[-1] is 0 because the boundary at max depth is 0)
        b[0] += (previous_time_lower_diagonal *
                 solution_grid[time_step - 1, 0] +
                 -current_time_lower_diagonal *
                 solution_grid[time_step, 0])
        # TESTING
        # b[0] = (-2 * previous_depth_factor * solution_grid[time_step, 0] +
        #         (8 * step_factor - current_depth_factor) *
        #         solution_grid[time_step - 1, 1] +
        #         (2 * step_factor - next_depth_factor) *
        #         solution_grid[time_step - 1, 2])
        # ###################

        # Calculate the current time step's solution
        # and insert it in between boundaries in solution matrix
        solution_grid[time_step, 1:] = crout_factorization(
            current_time_matrix, b)

    return solution_grid


def solomon(max_depth, max_time, depth_steps, time_steps,
            effective_diffusion_constant, effective_velocity_constant,
            surface_tracer_concentrations):
    """
    Uses the Crank Nicolson method to numerically approximate the PDE,
        du/dt = ED (d2u/dz2) - EV (du/dz)
    where the `d`s are partials, ED is the effective diffusion constant,
    EV is the effective velocity constant, z is depth, t is time, and u
    is the gas concentration.

    The surface tracer concentration is the bound for u(t, 0) where
    u(t, z). This must have a point for each step in the simulation and
    one for corner bound at u(0, 0).

    The boundaries along u(0, z) and u(t, max_depth + one_step) are
    assumed to be u = 0.

    Method
    ------
    The PDE is approximated using the central difference approximation
    for the second and first order partials of u with respect to z. Then
    using the Crank Nicolson method by using a weighted forward

    These approximations give a linear system with a tridiagonal matrix
    that is solved using the Crout Factorization algorithm.
    """
    assert len(surface_tracer_concentrations) == time_steps + 1, \
        'Surface boundary has length {}, not the expected {}'.format(
            len(surface_tracer_concentrations), time_steps + 1)
    delta_time = float(max_time) / time_steps
    delta_depth = float(max_depth) / depth_steps

    # The factors in iterative equation provided by the finite difference
    # approximations
    step_factor = np.power(delta_depth, 2) / delta_time
    # Note: These factors are for the current time step side of the equation
    # so the will be negated for calculations of the previous time step
    previous_depth_factor = (-6 * effective_diffusion_constant - 3 *
                             effective_velocity_constant * delta_depth)
    current_depth_factor = 12 * effective_diffusion_constant
    next_depth_factor = (-6 * effective_diffusion_constant + 3 *
                         effective_velocity_constant * delta_depth)

    # Construct the tridiagonal matrix for the current time step side
    current_time_upper_diagonal = 2 * step_factor + next_depth_factor
    current_time_middle_diagonal = 8 * step_factor + current_depth_factor
    current_time_lower_diagonal = 2 * step_factor + previous_depth_factor
    current_time_matrix = TridiagonalMatrix(
        depth_steps, [current_time_upper_diagonal] * (depth_steps - 1),
        [current_time_middle_diagonal] * (depth_steps),
        [current_time_lower_diagonal] * (depth_steps - 1))

    # Modify boundary factors
    # Despite appearing to do this Solomon's solver overwrites it
    # current_time_matrix.set(0, 0, 4 * step_factor - previous_depth_factor)
    current_time_matrix.set(
            current_time_matrix.size() - 1, current_time_matrix.size() - 1,
             4 * step_factor - previous_depth_factor)

    # Construct the tridiagonal matrix for the previous time step side
    previous_time_upper_diagonal = 2 * step_factor - next_depth_factor
    previous_time_middle_diagonal = 8 * step_factor - current_depth_factor
    previous_time_lower_diagonal = 2 * step_factor - previous_depth_factor
    previous_time_matrix = \
        (np.diag([previous_time_upper_diagonal] * (depth_steps - 1), 1) +
         np.diag([previous_time_middle_diagonal] * (depth_steps), 0) +
         np.diag([previous_time_lower_diagonal] * (depth_steps - 1), -1))

    # Modify boundary factors
    # Despite appearing to do this Solomon's solver overwrites it
    # previous_time_matrix[0, 0] = 4 * step_factor + previous_depth_factor
    previous_time_matrix[-1, -1] = 4 * step_factor + previous_depth_factor

    # Initialize the solution grid matrix. This matrix will contain all the
    # solution values, including the boundaries.
    solution_grid = np.empty((time_steps + 1, depth_steps + 1))
    # Fill with boundaries
    solution_grid[0,:] = 0
    solution_grid[:,0] = surface_tracer_concentrations
    for time_step in range(1, time_steps + 1):
        # Slice the boundary conditions from previous time step solution
        b = previous_time_matrix.dot(solution_grid[time_step - 1, 1:])

        # Add the boundary offset to b
        # (the offset for b[-1] is 0 because the boundary at max depth is 0)
        # b[0] += (previous_time_lower_diagonal *
        #          solution_grid[time_step - 1, 0] +
        #          current_time_lower_diagonal *
        #          solution_grid[time_step, 0])
        b[0] = (-2 * previous_depth_factor * solution_grid[time_step, 0] +
                (8 * step_factor - current_depth_factor) *
                solution_grid[time_step - 1, 1] +
                (2 * step_factor - next_depth_factor) *
                solution_grid[time_step - 1, 2])

        # Calculate the current time step's solution
        # and insert it in between boundaries in solution matrix
        solution_grid[time_step, 1:] = crout_factorization(
            current_time_matrix, b)

    return solution_grid


def thomas_algorithim(a, b, c, d):
    n = len(d)

    for k in range(1, n):
        m = a[k - 1] / b[k - 1]
        b[k] = b[k] - m * c[k - 1]
        d[k] = d[k] - m * d[k -1]

    x = np.empty(n)
    x[n - 1] = d[n - 1] / b[n - 1]

    for k in range(n - 2, -1, -1):
        x[k] = (d[k] - c[k] * x[k + 1]) / b[k]

    return x

def thomas_algorithim_single(a, b, c, d):
    n = len(d)

    b_prime = np.empty(n)
    b_prime[0] = b
    for k in range(1, n):
        m = a / b_prime[k - 1]
        b_prime[k] = b - m * c
        d[k] = d[k] - m * d[k -1]

    x = np.empty(n)
    x[n - 1] = d[n - 1] / b_prime[n - 1]

    for k in range(n - 2, -1, -1):
        x[k] = (d[k] - c * x[k + 1]) / b_prime[k]

    return x


def wieghted_time_single(max_depth, max_time, depth_steps, time_steps,
            effective_diffusion_constant, effective_velocity_constant,
            surface_tracer_concentrations):
    """
    Uses the Crank Nicolson method to numerically approximate the PDE,
        du/dt = ED (d2u/dz2) - EV (du/dz)
    where the `d`s are partials, ED is the effective diffusion constant,
    EV is the effective velocity constant, z is depth, t is time, and u
    is the gas concentration.

    The surface tracer concentration is the bound for u(t, 0) where
    u(t, z). This must have a point for each step in the simulation and
    one for corner bound at u(0, 0).

    The boundaries along u(0, z) and u(t, max_depth + one_step) are
    assumed to be u = 0.

    Method
    ------
    The PDE is approximated using the central difference approximation
    for the second and first order partials of u with respect to z. Then
    using the Crank Nicolson method by using a weighted forward

    These approximations give a linear system with a tridiagonal matrix
    that is solved using the Crout Factorization algorithm.
    """
    assert len(surface_tracer_concentrations) == time_steps + 1, \
        'Surface boundary has length {}, not the expected {}'.format(
            len(surface_tracer_concentrations), time_steps + 1)
    delta_time = float(max_time) / time_steps
    delta_depth = float(max_depth) / depth_steps

    # The factors in iterative equation provided by the finite difference
    # approximations
    step_factor = np.power(delta_depth, 2) / delta_time
    # Note: These factors are for the current time step side of the equation
    # so the will be negated for calculations of the previous time step
    previous_depth_factor = (-6 * effective_diffusion_constant - 3 *
                             effective_velocity_constant * delta_depth)
    current_depth_factor = 12 * effective_diffusion_constant
    next_depth_factor = (-6 * effective_diffusion_constant + 3 *
                         effective_velocity_constant * delta_depth)

    # Construct the tridiagonal matrix for the current time step side
    current_time_upper_diagonal = 2 * step_factor + next_depth_factor
    current_time_middle_diagonal = 8 * step_factor + current_depth_factor
    current_time_lower_diagonal = 2 * step_factor + previous_depth_factor

    # Construct the tridiagonal matrix for the previous time step side
    previous_time_upper_diagonal = 2 * step_factor - next_depth_factor
    previous_time_middle_diagonal = 8 * step_factor - current_depth_factor
    previous_time_lower_diagonal = 2 * step_factor - previous_depth_factor

    # Initialize the solution grid matrix. This matrix will contain all the
    # solution values, including the boundaries.
    solution_grid = np.empty((time_steps + 1, depth_steps + 1))
    # Fill with boundaries
    solution_grid[0,:] = 0
    solution_grid[:,0] = surface_tracer_concentrations
    previous_time_vector = np.empty(depth_steps)
    for time_step in range(1, time_steps + 1):
        previous_time_vector[0] = \
            (previous_time_middle_diagonal * solution_grid[time_step - 1, 1] +
             previous_time_upper_diagonal * solution_grid[time_step - 1, 2] +
             previous_time_lower_diagonal * solution_grid[time_step - 1, 0] -
             current_time_lower_diagonal * solution_grid[time_step, 0])
        for i in range(1, len(previous_time_vector) - 1):
            previous_time_vector[i] = \
                (previous_time_lower_diagonal * solution_grid[time_step - 1, i] +
                 previous_time_middle_diagonal * solution_grid[time_step - 1, i + 1] +
                 previous_time_upper_diagonal * solution_grid[time_step - 1, i + 2])
        previous_time_vector[len(previous_time_vector) - 1] = \
            (previous_time_lower_diagonal * solution_grid[time_step - 1, -2] +
             previous_time_middle_diagonal * solution_grid[time_step - 1, -1])

        # Calculate the current time step's solution
        # and insert it in between boundaries in solution matrix
        solution_grid[time_step, 1:] = thomas_algorithim_single(
            current_time_lower_diagonal, current_time_middle_diagonal,
            current_time_upper_diagonal, previous_time_vector)

    return solution_grid


def generate_solution_grid_csv(solver):
    d_star = 9.88011475
    q_star = 0.3
    theta_star = 0.2255
    effective_diffusion_constant = d_star / theta_star
    effective_velocity_constant = q_star / theta_star
    cfc_11_concentrations = np.genfromtxt(
        'cfc-11_atmospheric_concentrations.csv', dtype=float, delimiter=',',
        names=True)

    solution_grid = solver(
        200, 74, 1000, 148, effective_diffusion_constant,
        effective_velocity_constant, cfc_11_concentrations['concentration'])
    np.savetxt('solution_grid.csv', solution_grid, delimiter=',')

def plot_solution_grid_csv(csv):
    cfc_11_concentrations = np.genfromtxt(
        'cfc-11_atmospheric_concentrations.csv', dtype=float, delimiter=',',
        names=True)
    solution_grid = np.genfromtxt(
        csv, dtype=float, delimiter=',')

    for depth in range(0, 35, 5):
        plt.plot(cfc_11_concentrations['year'], solution_grid[:,depth * 5],
                 label=depth)

    plt.legend(loc='best', fancybox=True)
    plt.show()

def compare_solvers(first_solver, second_solver):
    d_star = 9.88011475
    q_star = 0.3
    theta_star = 0.2255
    effective_diffusion_constant = d_star / theta_star
    effective_velocity_constant = q_star / theta_star
    cfc_11_concentrations = np.genfromtxt(
        'cfc-11_atmospheric_concentrations.csv', dtype=float, delimiter=',',
        names=True)

    first_solution_grid = first_solver(
        200, 74, 1000, 148, effective_diffusion_constant,
        effective_velocity_constant, cfc_11_concentrations['concentration'])

    second_solution_grid = second_solver(
        200, 74, 1000, 148, effective_diffusion_constant,
        effective_velocity_constant, cfc_11_concentrations['concentration'])

    solution_grid_difference = first_solution_grid - second_solution_grid

    for depth in range(0, 35, 5):
        plt.plot(cfc_11_concentrations['year'], solution_grid_difference[:,depth * 5],
                 label=depth)

    plt.legend(loc='best', fancybox=True)
    plt.show()

def compare_csv_solution_grids(first_csv, second_csv):
    cfc_11_concentrations = np.genfromtxt(
        'cfc-11_atmospheric_concentrations.csv', dtype=float, delimiter=',',
        names=True)

    first_solution_grid = np.genfromtxt(first_csv, dtype=float, delimiter=',')
    second_solution_grid = np.genfromtxt(second_csv, dtype=float, delimiter=',')

    solution_grid_difference = (first_solution_grid - second_solution_grid) / second_solution_grid

    for depth in range(0, 35, 5):
        plt.plot(cfc_11_concentrations['year'], solution_grid_difference[:,depth * 5],
                 label=depth)

    plt.legend(loc='best', fancybox=True)
    plt.show()


if __name__ == '__main__':
    # compare_solvers(wieghted_time_single, wieghted_time)
    # generate_solution_grid_csv(wieghted_time_single)
    # compare_csv_solution_grids('./fi_solution_grid.csv', './solomon_fully_implicit.csv')
    # compare_csv_solution_grids('./cn_solution_grid.csv', './solomon_crank_nicolson.csv')
    compare_csv_solution_grids('./fi_solution_grid.csv', './cn_solution_grid.csv')
    # compare_csv_solution_grids('./solomon_fully_implicit.csv', './solomon_crank_nicolson.csv')
    # plot_solution_grid_csv('fi_solution_grid.csv')
