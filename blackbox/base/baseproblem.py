from abc import ABC, abstractmethod

class BaseProblem(ABC):
    """
        Abstract class for all the high level classes.
    """

    @abstractmethod
    def evaluate(self):
        """
            Method to perform a single run.
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
