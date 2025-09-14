from blackbox.airfoil_problem import AirfoilCST, AirfoilOptions
from pyDOE3 import lhs
import numpy as np

airfoil = AirfoilCST(options=AirfoilOptions())

airfoil.add_design_varaible("alpha", lower_bound=1.5, upper_bound=2.5)

airfoil.add_design_varaible("upper", lower_bound=np.array([-0.1]*6), upper_bound=np.array([0.1]*6))

x = lhs(airfoil.lower_bound.shape[0], 10, criterion="cm")

airfoil(x)
