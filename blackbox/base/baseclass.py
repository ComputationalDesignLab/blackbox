from abc import ABC, abstractmethod

class BaseClass(ABC):
    """
        Abstract class for all the high level classes.
    """

    def _error(self, message):
        """
            Method for printing errors in nice manner.
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

    def _getDefaultOptions(self, defaultOptions):
        """
            Setting up the initial values of options.

            Parameters
            ----------
            defaultOptions : dict
                Default options for the class.
        """

        for key in vars(defaultOptions):
            value = getattr(defaultOptions, key)
            self.options[key] = value

    def _setOptions(self, options):
        """
            Method for assigning user provided options.

            Parameters
            ----------
            options : dict
                User provided options.
        """

        for key in options.keys():
            # update strategy is different for a dictionary
            # and other type of variables
            if isinstance(options[key], dict):
                # if the option is already present in default dictionary, 
                # then add the user provided key-value pairs to the default dictionary.
                if key in self.options.keys():
                    self.options[key].update(options[key]) 
                else:
                    self.options[key] = options[key]
            else:
                self.options[key] = options[key]

    @abstractmethod
    def generateSamples(self):
        """
            Method to generate samples and save the data for further use.
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
