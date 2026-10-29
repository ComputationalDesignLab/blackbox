====================
Analytical Problems
====================

Blackbox is primarily designed for providing an interface for high-fidelity computationally expensively problems. However, for the sake of demonstration and completeness, Blackbox contains some analytical problems that are commonly used for benchmarking optimization algorithms and surrogate models. Currently, following problems are available for use:

* `Ackley function <https://www.sfu.ca/~ssurjano/ackley.html>`_
* `Levy function <https://www.sfu.ca/~ssurjano/levy.html>`_
* `Rastrigin function <https://www.sfu.ca/~ssurjano/rastr.html>`_
* `Hartmann 6D function <https://www.sfu.ca/~ssurjano/hart6.html>`_
* `Hartmann 4D function <https://www.sfu.ca/~ssurjano/hart4.html>`_
* `Hartmann 3D function <https://www.sfu.ca/~ssurjano/hart3.html>`_
* `Branin function <https://www.sfu.ca/~ssurjano/branin.html>`_
* Modified branin function

Basic usage
------------

To use an analytical problem provided in Blackbox, you can follow steps outlined below:

1. Initialize the problem class
2. Create a sampling plan (design of experiments)
3. Evaluate the sampling plan using the initialized problem class

The following code snippet shows an example of how to work with an analytical problem from Blackbox:

.. literalinclude:: ../../../examples/analytical_problems/ackley.py
    :start-after: # rst EXAMPLE start
    :end-before: # rst EXAMPLE end

Refer :doc:`api details <../api/analytical_problems>` for more information about how to use the analytical problems.