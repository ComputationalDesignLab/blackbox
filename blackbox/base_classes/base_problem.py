from abc import abstractmethod, ABC
import numpy as np

class BaseProblem(ABC):
    """
        An abstract base problem class which is inherited by all other problem classes.
    """

    @abstractmethod
    def __call__(self):
        
        pass

    def check_input(self, x: np.ndarray):
        """
            Method to validate the input before evaluation

            Parameters
            ----------
            x: np.ndarray
                Input sample(s) to be evaluated, can be either a 1D or 2D numpy array
        """

        x = np.atleast_2d(x)

        try:
            assert x.shape[1] == self.dim, "Input dimension must match the problem dimension"
            assert np.all(x >= self.lb) and np.all(x <= self.ub), "Input values must be within the bounds defined by lb and ub"
            
        except AssertionError as e:
            self._error(str(e))

    @staticmethod
    def _error(message: str, type=0) -> None:
        """
            Method for printing errors in nice format.

            Parameters
            ----------
            message: str
                Message to be displayed.

            type: int
                Type of the message - 0 for error and 1 for warning, default is 0.
        """

        # Initial message - total len is 80 characters
        msg = "\n+" + "-" * 78 + "+" + "\n" + "| Blackbox {}: ".format("Error" if type == 0 else "Warning")

        # Initial number of characters
        i = 16

        for word in message.split():
            if len(word) + i + 1 > 76:  # Finish line and start new one
                msg += " " * (76 - i) + " |\n| " + word + " " # Adding space and word in new line
                i = len(word) + 1 # Setting i value for new line
            else:
                msg += word + " " # Adding the word with a space
                i += len(word) + 1 # Increase the number of characters
        msg += " " * (76 - i) + " |\n" + "+" + "-" * 78 + "+" + "\n" # Adding last line
 
        print(msg, flush=True)

        if type == 0:
            exit()
