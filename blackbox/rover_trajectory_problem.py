# This module contains functions and classes used for defining rover trajectory
# and calculating the cost of the trajectory. The code in this module is based on
# the original code: https://github.com/zi-w/Ensemble-Bayesian-Optimization/tree/master/test_functions

from .base_classes.base_problem import BaseProblem
import numpy as np
import scipy.interpolate as si
import matplotlib.pyplot as plt

class Rover(BaseProblem):

    def __init__(
        self, 
        dim: int = 60, 
        start: np.ndarray = np.array([0.05, 0.05]), 
        goal: np.ndarray =np.array([0.95, 0.95]),
        force_start: bool = False,
        force_goal: bool = False
    ):
        """
            Main class for the user to create a domain with obstacles 
            and evaluate the cost of a trajectory for a given input

            Parameters
            -----------
            dim: int
                dimension of the problem, determines the number of control points for the B-spline trajectory

            start: np.ndarray
                1D numpy array representing (x,y) coordinates of the starting point, default=(0.05,0.05)

            goal: np.ndarray
                1D numpy array representing (x,y) coordinates of the goal point, default=(0.95,0.95)

            force_start: bool
                flag to fix the starting point, default=False

            force_goal: bool
                flag to fix the goal point, default=False
        """

        assert isinstance(dim, int), "dim should of int type"
        assert isinstance(start, np.ndarray), "start should be a numpy array"
        assert isinstance(goal, np.ndarray), "goal should be a numpy array"
        assert isinstance(force_start, bool), "force_start should be a bool"
        assert isinstance(force_goal, bool), "force_goal should be a bool"

        # some parameters
        self.dim = dim
        self.lb = np.zeros(dim)
        self.ub = np.ones(dim)
        env_dim = 2 # dimension of the domain

        # n_points = number of control points for the B-spline trajectory
        num_control_points = int(dim/env_dim)

        # create cost function
        cost_fn = create_cost()

        # create trajectory
        traj = PointBSpline(dim=env_dim, num_control_points=num_control_points)

        # create domain
        self.domain = RoverDomain(
            cost_fn,
            start=start,
            goal=goal,
            traj=traj,
            start_miss_cost=l1cost,
            goal_miss_cost=l1cost,
            force_start=force_start,
            force_goal=force_goal
        )

        # range of states i.e. range of each (x,y) coordinate - NOT sure why this is needed
        s_range = np.array([[-0.1, -0.1], [1.1, 1.1]])

        # lb and ub for all the design variable
        raw_x_range = np.repeat(s_range, num_control_points, axis=1)

        # minimum value of f - from the original code
        f_min = 5.0

        # setup the function handle to call for evaluation
        # not sure why normalizing the input to s_range - borrowing from original code
        self._func = NormalizedInputFn(ConstantOffsetFn(self.domain, f_min), raw_x_range)

    def __call__(self, x: np.ndarray) -> np.ndarray:
        """
            Evaluate the cost of the trajectory for a given x

            Parameters
            ----------
            x : np.ndarray
                input to be evaluated, can be a single sample or multiple samples, (dim,) or (n_samples,dim)

            Returns
            -------
            cost : np.ndarray
                cost of the trajectory for the given input, (1,) or (n_samples,1)
        """

        assert isinstance(x, np.ndarray), f"x should be a numpy array, but got {type(x)}"

        if x.ndim == 1:
            x = x.reshape(1,-1)

        y = np.zeros((x.shape[0],1))

        for i in range(x.shape[0]):
            y[i,0] = self._func(x[i,:])

        if x.ndim == 1:
            y = y.reshape(-1,)

        return y


class NormalizedInputFn:

    def __init__(self, fn_instance, x_range):
        self.fn_instance = fn_instance
        self.x_range = x_range

    def __call__(self, x):
        return self.fn_instance(self.project_input(x))

    def project_input(self, x):
        return x * (self.x_range[1] - self.x_range[0]) + self.x_range[0]

    def inv_project_input(self, x):
        return (x - self.x_range[0]) / (self.x_range[1] - self.x_range[0])

    def get_range(self):
        return np.array([np.zeros(self.x_range[0].shape[0]), np.ones(self.x_range[0].shape[0])])
    

class ConstantOffsetFn:

    def __init__(self, fn_instance, offset):
        """
            Class to add a constant offset to a given function's output

            Parameters
            -----------
            fn_instance: function handle

            offset: float
                offset to be added
        """
        
        self.fn_instance = fn_instance
        self.offset = offset

    def __call__(self, x):
        """
            Method to return the function output with an offset
        """

        return self.fn_instance(x) + self.offset
    

################################################
############## Trajectory classes ##############
################################################

class PointBSpline:

    def __init__(self, dim: int, num_control_points: int):
        """
            Class to create a smooth B-spline trajectory based on the given points

            Parameters
            ----------

            dim: int
                number of dimensions of the domain (env), mostly 2D

            num_control_points: int
                number of control points used for fitting B-spline trajectory
        """
        
        self.dim_domain = dim
        self.num_control_points = num_control_points
        self.tck = None # variable to represent knot vector, coefficients, and degree of the b-spline

    def fit_bspline(self, control_points: np.ndarray, start: np.ndarray = None, goal: np.ndarray = None):
        """
            Method to fit a parameterized b-spline curve

            Parameters
            ----------
            control_points: np.ndarray
                1D numpy array containing (x,y) coordinates of the control points

            start: np.ndarray
                1D numpy array containing (x,y) for the starting position, default=None

            goal: np.ndarray
                1D numpy array containing (x,y) for the goal position, default=None
        """

        control_points = control_points.reshape(-1,self.dim_domain).T

        if start is not None:
            control_points = np.hstack(( start.reshape(-1,1), control_points ))

        if goal is not None:
            control_points = np.hstack(( control_points, goal.reshape(-1,1) ))

        # fit the b-spline, k = degree
        self.tck, _ = si.splprep(control_points, k=3)

        # Setting coefficients to start and goal value
        if start is not None:
            for a, sv in zip(self.tck[1], start):
                a[0] = sv

        if goal is not None:
            for a, gv in zip(self.tck[1], goal):
                a[-1] = gv

    def get_traj_points(self, u: np.ndarray) -> np.ndarray:
        """
            Method to evaluate the parameterized b-spline curve at given
            parameter values

            Parameters
            ----------
            u: np.ndarray
                1D numpy array consisting of parameter values for the b-spline curve, shape (n,)
                where n is the number of parameter values

            Returns
            --------
            traj_points: np.ndarray
               a numpy array representing points along the b-spline curve, shape (n,d)
        """

        assert self.tck is not None, "Fit the B-spline before evaluating it"
        assert u.ndim == 1

        traj_points = np.vstack(si.splev(u, self.tck)).T

        return traj_points
    
################################################
############# Rover Domain class ###############
################################################

class RoverDomain():

    def __init__(self,
        cost_fn,
        start: np.ndarray,
        goal: np.ndarray,
        traj,
        start_miss_cost=None,
        goal_miss_cost=None,
        force_start: bool = True,
        force_goal: bool = True,
    ):
        """
            Class for creating rover domain

            Parameters
            ----------
            cost_fn: 
                function handle to compute cost for a given set of points that define a trajectory

            start: np.ndarray
                starting position (x,y) of the rover, shape (2,)

            end: np.ndarray
                final (end) position (x,y) of the rover, shape (2,)

            traj:
                a parameterized trajectory object that can create a smooth b-spline trajectory
                from the given control points

            start_miss_cost:
                function handle to compute cost for missing starting position, default=None

            goal_miss_cost:
                function handle to compute cost for missing final position, default=None
                
            force_start: bool
                flag to compute cost for missing the start position, default=True

            force_goal: bool
                flag to compute cost for missing the goal position, default=True
        """

        assert start.ndim == 1
        assert goal.ndim == 1

        self.cost_fn = cost_fn
        self.start = start
        self.goal = goal
        self.dim = start.shape[0]
        self.traj = traj
        self.start_miss_cost = start_miss_cost
        self.goal_miss_cost = goal_miss_cost
        self.force_start = force_start
        self.force_goal = force_goal

        # assign cost function for missing start or end position, if not provided by user
        if self.start_miss_cost is None:
            self.start_miss_cost = simple_rbf
        if self.goal_miss_cost is None:
            self.goal_miss_cost = simple_rbf

    def __call__(self, control_points: np.ndarray, n_samples: int = 1000) -> float:
        """
            Method to estimate the cost for a given control points i.e.
            design variable vector x

            Parameters
            ----------
            control_points: np.ndarray
                1D numpy array representing coordinates of the 
                control points for fitting b-spline, shape (n,) where
                n = number of control points X 2. This is essential the design
                variable vector for which cost needs to be computed

            n_samples: int
                number of points which will be used to query the 
                b-spline after fitting, default = 1000

            Returns
            -------
            total_cost: float
                total cost of the trajectory
        """

        # fit b-spline for given control points
        # Not sure why the noise is added to the control point
        self.traj.fit_bspline(control_points + np.random.normal(loc=0,scale=1e-4,size=control_points.shape),
                             self.start if self.force_start else None,
                             self.goal if self.force_goal else None)
        
        # query points for the parameterized b-spline
        u = np.linspace(0.0, 1.0, n_samples, endpoint=True)

        # points on the trajectory
        traj_points = self.traj.get_traj_points(u)

        costs = self.cost_fn(traj_points).reshape(-1,) # cost at each point - shape (n,)

        # estimate cost for the entire trajectory using trapezoidal integration
        # Note: first and last point is not included in the integration
        avg_cost = 0.5 * (costs[1:] + costs[:-1])
        delx = np.linalg.norm(traj_points[1:] - traj_points[:-1], axis=1)
        total_cost = np.sum(delx * avg_cost)

        # total cost if the start and goal are part of objective
        if not self.force_start:
            total_cost += self.start_miss_cost(traj_points[0], self.start)
        if not self.force_goal:
            total_cost += self.goal_miss_cost(traj_points[-1], self.goal)

        return total_cost
    

class AABoxes():
    
    def __init__(self, low: np.ndarray, high: np.ndarray):
        """
            Class to create axis aligned bounding boxes for the obstacles,
            typically used to model rectangular obstacles in 2D or 3D spaces

            Parameters
            ----------
            low : np.ndarray
                lower bound of the bounding box, should be a 2D array of shape (n, dim)
                where n is the number of boxes and dim is the dimension of the space

            high : np.ndarray
                upper bound of the bounding box, should be a 2D array of shape (n, dim)
                where n is the number of boxes and dim is the dimension of the space
        """

        self.low = low
        self.high = high

    def contains(self, X: np.ndarray) -> np.ndarray:
        """
            Check if a set of points X is inside the bounding box

            Parameters
            ----------
            X : np.ndarray
                2D numpy array containing (x,y) coordinates

            Returns
            -------
            bool : np.ndarray
                2D boolean numpy array indicating if the point is inside the bounding box
        """

        assert X.ndim == 2, "X should be 2D numpy array"
    
        # convert the arrays to higher dimension so that numpy
        # arrays of different sizes can be compared
        lX = self.low.T[None, :, :] <= X[:, :, None]
        hX = self.high.T[None, :, :] > X[:, :, None]

        return (lX.all(axis=1) & hX.all(axis=1))
    

class NegGeom():

    def __init__(self, geom):
        """
            Class to create a negation of the given geometry. The negation
            is used to transform the occupied regions into free space
            and vice versa, which is particularly useful in scenarios 
            where defining the complement of a geometric region is more 
            straightforward or efficient.

            Parameters
            ----------
            geom : object
                geometry to be negated, should be a class that implements the
                `contains` method to check if a point is inside the geometry
        """

        self.geom = geom

    def contains(self, X: np.ndarray) -> np.ndarray:
        """
            Check if a set of points X is inside the negated geometry. 

            Parameters
            ----------
            x : np.ndarray
                2D numpy array containing (x,y) coordinates

            Returns
            -------
            bool : np.ndarray
                2D boolean numpy array indicating if the point is inside the bounding box
        """

        return ~self.geom.contains(X)
    

class UnionGeom():

    def __init__(self, geoms: list):
        """
            Class to create a union of the given geometries. The union
            is used to combine multiple geometric regions into a single
            region, which is useful for defining complex obstacles or
            free space.

            Parameters
            ----------
            geoms : list
                list of geometries to be combined, should be a list of classes 
                that implements the `contains` method to calculate the cost of the trajectory
        """

        self.geoms = geoms

    def contains(self, X: np.ndarray) -> np.ndarray:
        """
            Check if the point X is inside the union of the geometries.

            Parameters
            ----------
            X : np.ndarray
                2D numpy array containing (x,y) coordinates

            Returns
            -------
            bool : np.ndarray
                2D boolean array indicating if the point is inside the union of the geometries
        """

        return np.any(np.hstack([geom.contains(X) for geom in self.geoms]), axis=1, keepdims=True)
    

################################################
################ Cost functions ################
################################################

class ConstObstacleCost():

    def __init__(self, geometry, obstacle_cost: float):
        """
            class to create a cost function for the obstacles

            Parameters
            ----------
            geometry : object
                geometry of the obstacles, should be a class that implements the
                __call__ method to calculate the cost of the trajectory

            obstacle_cost : float
                constant obstacle cost if a given point is inside the obstacles
        """
        
        self.geom = geometry
        self.obstacle_cost = obstacle_cost

    def __call__(self, x: np.ndarray) -> np.ndarray:
        """
            method to compute obstacle for each point

            Parameters
            ----------
            x: np.ndarray
                numpy array containing (x,y) point(s)

            Returns
            -------
            2D numpy array containing constant obstacle cost for each point
        """

        return self.geom.contains(x) * self.obstacle_cost
    

class ConstCost():

    def __init__(self, constant_cost: float):
        """
            class to create a constant cost function

            Parameters
            ----------
            cost : float
                constant cost for each point in the trajectory
        """
        
        self.constant_cost = constant_cost

    def __call__(self, x: np.ndarray) -> np.ndarray:
        """
            Method to assign a constant cost for each point

            Parameters
            ----------
            x: np.ndarray
                numpy array containing (x,y) point(s)

            Returns
            -------
            2D numpy array containing constant cost for each point
        """

        if x.ndim == 1:
            x = x.reshape(1, -1)

        return self.constant_cost * np.ones((x.shape[0], 1))
    

class AdditiveCosts():
    
    def __init__(self, fns: list):
        """
            Class to create a cost function that is the sum of the given cost functions.

            Parameters
            ----------
            fns : list
                list of cost function handles to be added
        """
        
        self.fns = fns

    def __call__(self, x: np.ndarray) -> np.ndarray:
        """
            Method to compute additive cost for given x

            Parameters
            ----------
            x: np.ndarray
                numpy array containing (x,y) point(s)

            Returns
            -------
            2D numpy array containing constant cost for each point
        """

        return np.sum([fn(x) for fn in self.fns], axis=0)

def create_cost():
    """
        function to create a cost function, the location
        of the obstacles is fixed to ensure consistency

        Returns
        -------
        cost_fn : function
            function handle for calculating the cost of the trajectory
    """

    # center of the obstacles
    c = np.array([[0.43143755, 0.20876147],
                  [0.38485367, 0.39183579],
                  [0.02985961, 0.22328303],
                  [0.7803707, 0.3447003],
                  [0.93685657, 0.56297285],
                  [0.04194252, 0.23598362],
                  [0.28049582, 0.40984475],
                  [0.6756053, 0.70939481],
                  [0.01926493, 0.86972335],
                  [0.5993437, 0.63347932],
                  [0.57807619, 0.40180792],
                  [0.56824287, 0.75486851],
                  [0.35403502, 0.38591056],
                  [0.72492026, 0.59969313],
                  [0.27618746, 0.64322757],
                  [0.54029566, 0.25492943],
                  [0.30903526, 0.60166842],
                  [0.2913432, 0.29636879],
                  [0.78512072, 0.62340245],
                  [0.29592116, 0.08400595],
                  [0.87548394, 0.04877622],
                  [0.21714791, 0.9607346],
                  [0.92624074, 0.53441687],
                  [0.53639253, 0.45127928],
                  [0.99892031, 0.79537837],
                  [0.84621631, 0.41891986],
                  [0.39432819, 0.06768617],
                  [0.92365693, 0.72217512],
                  [0.95520914, 0.73956575],
                  [0.820383, 0.53880139],
                  [0.22378049, 0.9971974],
                  [0.34023233, 0.91014706],
                  [0.64960636, 0.35661133],
                  [0.29976464, 0.33578931],
                  [0.43202238, 0.11563227],
                  [0.66764947, 0.52086962],
                  [0.45431078, 0.94582745],
                  [0.12819915, 0.33555344],
                  [0.19287232, 0.8112075],
                  [0.61214791, 0.71940626],
                  [0.4522542, 0.47352186],
                  [0.95623345, 0.74174186],
                  [0.17340293, 0.89136853],
                  [0.04600255, 0.53040724],
                  [0.42493468, 0.41006649],
                  [0.37631485, 0.88033853],
                  [0.66951947, 0.29905739],
                  [0.4151516, 0.77308712],
                  [0.55762991, 0.26400156],
                  [0.6280609, 0.53201974],
                  [0.92727447, 0.61054975],
                  [0.93206587, 0.42107549],
                  [0.63885574, 0.37540613],
                  [0.15303425, 0.57377797],
                  [0.8208471, 0.16566631],
                  [0.14889043, 0.35157346],
                  [0.71724622, 0.57110725],
                  [0.32866327, 0.8929578],
                  [0.74435871, 0.47464421],
                  [0.9252026, 0.21034329],
                  [0.57039306, 0.54356078],
                  [0.56611551, 0.02531317],
                  [0.84830056, 0.01180542],
                  [0.51282028, 0.73916524],
                  [0.58795481, 0.46527371],
                  [0.83259048, 0.98598188],
                  [0.00242488, 0.83734691],
                  [0.72505789, 0.04846931],
                  [0.07312971, 0.30147979],
                  [0.55250344, 0.23891255],
                  [0.51161315, 0.46466442],
                  [0.802125, 0.93440495],
                  [0.9157825, 0.32441602],
                  [0.44927665, 0.53380074],
                  [0.67708372, 0.67527231],
                  [0.81868924, 0.88356194],
                  [0.48228814, 0.88668497],
                  [0.39805433, 0.99341196],
                  [0.86671752, 0.79016975],
                  [0.01115417, 0.6924913],
                  [0.34272199, 0.89543756],
                  [0.40721675, 0.86164495],
                  [0.26317679, 0.37334193],
                  [0.74446787, 0.84782643],
                  [0.55560143, 0.46405104],
                  [0.73567977, 0.12776233],
                  [0.28080322, 0.26036748],
                  [0.17507419, 0.95540673],
                  [0.54233783, 0.1196808],
                  [0.76670967, 0.88396285],
                  [0.61297539, 0.79057776],
                  [0.9344029, 0.86252764],
                  [0.48746839, 0.74942784],
                  [0.18657635, 0.58127321],
                  [0.10377802, 0.71463978],
                  [0.7771771, 0.01463505],
                  [0.7635042, 0.45498358],
                  [0.83345861, 0.34749363],
                  [0.38273809, 0.51890558],
                  [0.33887574, 0.82842507],
                  [0.02073685, 0.41776737],
                  [0.68754547, 0.96430979],
                  [0.4704215, 0.92717361],
                  [0.72666234, 0.63241306],
                  [0.48494401, 0.72003268],
                  [0.52601215, 0.81641253],
                  [0.71426732, 0.47077212],
                  [0.00258906, 0.30377501],
                  [0.35495269, 0.98585155],
                  [0.65507544, 0.03458909],
                  [0.10550588, 0.62032937],
                  [0.60259145, 0.87110846],
                  [0.04959159, 0.535785]])
    
    # low and high values of the obstacles
    l = c - 0.025
    h = c + 0.025

    # boundary of the env
    r_box = np.array([[0.5, 0.5]]) # center of the env
    r_l = r_box - 0.5
    r_h = r_box + 0.5

    # create obstalces
    trees = AABoxes(l, h)

    # create the env
    r_box = NegGeom(AABoxes(r_l, r_h))

    # create obstacles
    obstacles = UnionGeom([trees, r_box])

    costs = [ConstObstacleCost(obstacles, obstacle_cost=20.), ConstCost(0.05)] # list of cost functions
    cost_fn = AdditiveCosts(costs) # add all the costs

    return cost_fn

def l1cost(x, point):
    """
        function to calculate the l1 cost if the x deviates from point

        Parameters
        ----------
        x: np.ndarray
            point to be evaluated
        point: np.ndarray
            point to be matched

        Returns
        -------
        a float representing L1 cost
    """

    return 10 * np.linalg.norm(x - point, 1)

def simple_rbf(x, point):
    """
        function to calculate the RBF-based cost if the x deviates from point

        Parameters
        ----------
        x: np.ndarray
            point to be evaluated
        point: np.ndarray
            point to be matched

        Returns
        -------
        a float representing L1 cost
    """
    
    return (1 - np.exp(-np.sum(((x - point) / 0.25) ** 2)))

def plot_traj(control_points, rover, ngrid_points=100, ntraj_points=1000, colormap='GnBu', draw_colorbar=True, show_figure=True):
    """
        function to plot the trajectory for given control points

        Parameters
        ----------
        control_points: np.ndarray
            a 1D numpy array representing the coordinates of the control points
        
        rover: Rover
            a Rover object

        ngrid_points: int
            number of grid points along one axis, default=100

        ntraj_points: int
            number of points along the trajectory, default=1000

        colormap: str
            a valid matplotlib colormap, default="GnBu"

        draw_colorbar: bool
            flag to have colorbar or not, default=True

        show_figure: bool
            flag to determine whether to show the plot or not, default=True
    """

    total_cost = rover(control_points)

    # get a grid of points over the state space
    s_range = np.array([[-0.1, -0.1], [1.1, 1.1]])
    points = [np.linspace(mi, ma, ngrid_points, endpoint=True) for mi, ma in zip(*s_range)]
    grid_points = np.meshgrid(*points)
    points = np.hstack([g.reshape((-1, 1)) for g in grid_points])

    # compute the cost at each point on the grid - not for trajectory
    costs = rover.domain.cost_fn(points)

    # get points on the current trajectory
    traj_points = rover.domain.traj.get_traj_points(np.linspace(0., 1.0, ntraj_points, endpoint=True))

    fs = 12 # fontsize

    fig, ax = plt.subplots(figsize=(8,6))

    # set title to be the total cost
    ax.set_title(f'Trajectory cost: {total_cost}')
    print(total_cost)

    # plot cost function
    cmesh = ax.pcolormesh(grid_points[0], grid_points[1], costs.reshape((ngrid_points, -1)), cmap=colormap)

    if draw_colorbar:
        cbar = plt.gcf().colorbar(cmesh)
        cbar.ax.set_ylabel("$c(x)$", fontsize=fs)

    # plot traj
    ax.plot(traj_points[:, 0], traj_points[:, 1], 'r')

    # plot start and goal
    ax.scatter([rover.domain.start[0], rover.domain.goal[0]], [rover.domain.start[1], rover.domain.goal[1]], marker="o", c="k")

    ax.tick_params(axis='both', labelsize=fs-2)
    plt.tight_layout()

    # show plot
    if show_figure:
        plt.show()

    # plt.close()

    return fig, ax
