=========
Blackbox
=========

Comples simulations are ubiquitous in various engineering and scientific domains. 

These simulations often involve solving complex set of equations, running time-consuming numerical methods, or executing multi-physics models.

This makes it harder to directly use these simulations for tasks such as optimization, uncertainty quantification, sensitivity analysis, etc.

To alleviate this challenge, surrogate models are often employed to approximate the behavior of these simulations.

These models are built using the data generated from running these simulations at given input points.

Typically, running these simulations requires significant efforts in terms of pre-processing, configuration, simulation execution, and post-processing.

The primary motivation behind **Blackbox** is to solve this issue by providing a consistent, easy-to-use API for evaluating such problems for a given input.

By abstracting away the complexities of problem/simulation setup, execution, and data handling, the package allows users to focus on higher-level tasks such as surrogate 
modeling, optimization, uncertainty quantification, etc., rather than on managing individual simulation workflows.

The user provides samples (or inputs) for evaluation to Blackbox, which internally handles the entire workflow and returns the desired outputs.

This is illustrated in the below figure as well.

.. image:: _static/bb.png
    :alt: blackbox motivation image

For example, if user wants to optimize the performance of an airfoil in unsteady flow using Bayesian optimization.



Blackbox is being actively developed to contain a wide range of problems from different domains.

These problems range from simple analytical functions to complex engineering simulation pipelines. 

Currently, following problems are implemented within Blackbox:

.. .. toctree::
..    :maxdepth: 3

..    install
..    airfoil/intro
..    wing/intro
..    hpc