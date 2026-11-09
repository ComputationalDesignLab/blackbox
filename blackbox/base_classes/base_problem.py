from abc import abstractmethod, ABC
import numpy as np

class BaseProblem(ABC):
    """
        An abstract base problem class which is inherited by all other problem classes.
    """

    @abstractmethod
    def __call__(self, x: np.ndarray) -> np.ndarray:
        pass
