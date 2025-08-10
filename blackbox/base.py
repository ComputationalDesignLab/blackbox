from abc import abstractmethod, ABC

class BaseFunction(ABC):
    """
        Base class providing a template for other model classes
    """

    @abstractmethod
    def evaluate(self, x):
        """
            Method for computing function values for given x.

            Parameters
            ----------
            x: numpy array
                numpy array containing  samples to be evaluated. x can be 
                a single sample or multiple samples
        """

        pass
