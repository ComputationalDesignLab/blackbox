from abc import ABC, abstractmethod

class BaseProblem(ABC):
    """
        Abstract class for all the high level classes.
    """

    @abstractmethod
    def evaluateSample(self):
        """
            Method to evaluate a single sample.
        """

        pass

    @abstractmethod
    def evaluateDOE(self):
        """
            Method to evaluate a DOE/Sampling plan.
        """

        pass

    @abstractmethod
    def addDV(self):
        """
            Method to add design variables.
        """

        pass

    @abstractmethod
    def removeDV(self):
        """
            Method to remove design variables.
        """

        pass
