# This module defines various function for generating sampling plan using different methods

from pyDOE3 import lhs
from scipy.stats.qmc import Sobol, Halton

def generate_lhs_samples(problem, n_samples, criterion="cm", iterations=100, random_state=None):

    samples = lhs(problem.in_dim, n_samples, criterion, iterations, random_state)

    samples = problem.lb + (problem.ub - problem.lb) * samples

    return samples

def generate_sobol_samples(problem, n_samples, scramble=False, random_state=None):

    pass

def generate_halton_samples(problem, n_samples, scramble=False, random_state=None):

    pass
