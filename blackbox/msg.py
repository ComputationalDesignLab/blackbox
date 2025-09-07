# This module contains methods for printing error in a formatted way for the blackbox package

def print_msg(message: str, type=0) -> None:
    """
        Method for printing blackbox messages in nice format.

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
