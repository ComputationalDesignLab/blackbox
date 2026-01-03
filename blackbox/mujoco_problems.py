# This module contains a class for defining and running a Mojoco environment
# Most of the code in this module is adapted from: 
# https://github.com/LamNgo1/cma-meta-algorithm/blob/master/test_functions/function_realworld_bo/functions_mujoco.py

import numpy as np
from typing import ClassVar, Tuple
from .base_classes.base_problem import BaseProblem

try:
    import gymnasium as gym
except ImportError:
    raise ImportError("gymnasium is not installed. Please install it to use Mujoco environments.")

class RunningStat():

    def __init__(self, shape=None):
        """
            Class for computing statistics for a stream of data,
            here it is used for computing the mean and standard deviation
            of the observation space of the Mujoco env. This mean and 
            standard deviation is used for normalizing the observation space
            before computing the next action

            http://www.johndcook.com/blog/standard_deviation/
        """

        self._n = 0  # Number of samples
        self._M = np.zeros(shape, dtype=np.float64) # mean
        self._S = np.zeros(shape, dtype=np.float64) # sum of squares

    def push(self, x: np.ndarray):
        """
            Method for adding a new value of x i.e new observation state
            to the running statistics

            Note: x should be a single observation state, not a batch of states
        """

        assert x.shape == self._M.shape, ("x.shape = {}, self.shape = {}".format(x.shape, self._M.shape))

        n1 = self._n
        self._n += 1

        if self._n == 1:
            self._M = x

        else:
            delta = x - self._M
            self._M += delta/self._n
            self._S += delta*delta*n1/self._n

    @property
    def n(self):
        return self._n

    @property
    def mean(self):
        return self._M

    @property
    def var(self):
        return self._S / (self._n - 1) if self._n > 1 else np.square(self._M)

    @property
    def std(self):
        return np.sqrt(self.var)
    
    @property
    def shape(self):
        return self._M.shape


class MujocoPolicyFunc(BaseProblem):

    # class variables
    SWIMMER_ENV: ClassVar[Tuple[str, float, float, int]] = ('Swimmer-v5', -1.0, 1.0, 3)
    HALF_CHEETAH_ENV: ClassVar[Tuple[str, float, float, int]] = ('HalfCheetah-v5', -1.0, 1.0, 3)
    HOPPER_ENV: ClassVar[Tuple[str, float, float, int]] = ('Hopper-v5', -1.0, 1.0, 3)
    WALKER_2D_ENV: ClassVar[Tuple[str, float, float, int]] = ('Walker2d-v5', -1.0, 1.0, 3)
    HUMANOID_ENV: ClassVar[Tuple[str, float, float, int]] = ('Humanoid-v5', -1.0, 1.0, 3)

    def __init__(self, env_name: str, lb: float, ub: float, num_rollouts: int):
        """
            Base class for all mujoco environments

            Parameters
            ----------
            env_name : str
                name of the mujoco environment
            lb : float
                lower bound of the input space
            ub : float
                upper bound of the input space
            num_rollouts : int
                number of times to run the entire episode for a given policy
        """

        self.env_name = env_name
        self.env = gym.make(env_name)
        self.env.reset(seed=2025) # imp to set this seed for consistency

        state_dim = self.env.observation_space.shape[0]
        action_dim = self.env.action_space.shape[0]
        self.dim = state_dim * action_dim
        self.policy_shape = (action_dim, state_dim)

        self.lb = np.full(self.dim, lb)
        self.ub = np.full(self.dim, ub)
        self.num_rollouts = num_rollouts # no. of times to run the entire episode for a given policy

    def __call__(self, x: np.ndarray) -> np.ndarray:
        """
            Method to compute total reward for a given policy x

            Parameters
            ----------
            x : np.ndarray
                input to be evaluated, can be a single sample or multiple samples, (dim,) or (n_samples,dim)

            Returns
            -------
            reward : np.ndarray
                reward of the trajectory for the given input, (1,) or (n_samples,1)
        """

        assert isinstance(x, np.ndarray), f"x should be a numpy array, but got {type(x)}"

        if x.ndim == 1:
            x = x.reshape(1,-1)

        assert np.all(x <= self.ub) and np.all(x >= self.lb)
        assert x.shape[1] == self.dim, f"x.shape[1] = {x.shape[1]}, expected {self.dim}"

        y = np.zeros((x.shape[0],1))

        for i in range(x.shape[0]):

            M = x[i,:].reshape(self.policy_shape) # policy matrix
            self.rs = RunningStat(self.env.observation_space.shape[0]) # initialize running stats class
            total_reward = 0.0

            for _ in range(self.num_rollouts):

                obs, _ = self.env.reset()

                while True:

                    self.rs.push(obs)

                    norm_obs = (obs - self.rs.mean) / (self.rs.std + 1e-6) # normalize observation

                    action = np.dot(M, norm_obs) # next action

                    obs, reward, terminated, truncated, _ = self.env.step(action) # compute next observation and reward

                    total_reward += reward

                    if terminated or truncated:
                        break

            y[i,0] = -total_reward / self.num_rollouts

        if x.ndim == 1:
            y = y.reshape(-1,)

        return y

class HalfCheetah(MujocoPolicyFunc):

    def __init__(self):
        super().__init__(*MujocoPolicyFunc.HALF_CHEETAH_ENV)