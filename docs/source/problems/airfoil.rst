=================
Airfoil Problem
=================

The ``Airfoil`` problem provides an interface for generating aerodynamic datasets for airfoils using high-fidelity computational fluid dynamics (CFD) workflow. Users can parameterize both the airfoil geometry and flow conditions and evaluate the resulting configurations without manually setting up the geometry, meshing, solver, and post-processing workflow. The ``Airfoil`` problem currently uses ADflow_ as the CFD solver and the Class-Shape Transformation (CST_) method for airfoil geometry parameterization. Additional CFD solvers, such as OpenFOAM and SU2, may be supported in future releases. The following sections describe the dependencies, computational workflow, and how to generate aerodynamic datasets with the ``Airfoil`` problem.

Dependencies
=============

The ``Airfoil`` problem requires several additional packages that must be installed before using this problem. The required packages and their corresponding versions are listed below.

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Package
     - Version
     - Comments
   * - `mpi4py`_
     - 3.1.6
     - Used for parallel processing during volume mesh generation and CFD solver execution
   * - `mdolab-baseclasses`_
     - >= 1.8.4
     - Provides ``AeroProblem`` for defining flow and operating conditions
   * - `pyHyp`_
     - >= 2.6.2
     - Used for generating the CFD volume mesh from the surface mesh
   * - `ADflow`_
     - 2.11.0
     - Used as the CFD solver for airfoil flow simulations
   * - `pyvista`_
     - `-`
     - Used for reading, processing, and writing surface and volume field data
   * - `h5py`_
     - `-`
     - Used for creating and storing field results in HDF5 format

.. _CST: https://doi.org/10.2514/1.29958
.. _mpi4py: https://mpi4py.readthedocs.io/en/3.1.6/
.. _mdolab-baseclasses: https://github.com/mdolab/baseclasses
.. _pyHyp: https://github.com/mdolab/pyhyp
.. _ADflow: https://github.com/mdolab/adflow
.. _pyvista: https://docs.pyvista.org/
.. _h5py: https://docs.h5py.org/

Make sure that all required packages are installed with compatible versions before using the Airfoil problem.

Workflow
=========

The ``Airfoil`` problem accepts a set of samples, where each sample contains parameter values. This set may contain a single sample or multiple samples. For each sample :math:`x_i` in the sample set, the ``Airfoil`` problem performs the sequence of operations described below:

1. **Airfoil shape parameterization**: The geometry parameters are used to define the airfoil shape through CST parameterization
2. **Surface mesh generation**: The resulting airfoil coordinates are used to generate the surface mesh
3. **Volume mesh generation**: A volume mesh is generated from the surface mesh using ``pyHyp``
4. **CFD analysis**: The volume mesh and the flow parameters are passed to the ``ADflow`` solver to compute the aerodynamic solution
5. **Output extraction**: Scalar and field quantities from the solver are extracted and stored in a folder

When multiple samples are provided, the workflow is evaluated for each sample, and the resulting outputs are stored separately for each sample. This allows the Airfoil problem to provide a consistent interface for generating datasets without requiring users to manually manage geometry processing, mesh generation, solver execution, or output handling.

The complete workflow is illustrated in the XDSM diagram below.

.. image:: ../_static/xdsm_airfoil.png
    :width: 100%
    :align: center
    :alt: airfoil problem xdsm
    :class: image-margin

The stacked rectangles representing the volume-mesh generation and CFD analysis indicate that these operations use multiple processors through MPI. The ``Airfoil`` problem therefore abstracts the geometry processing, mesh generation, solver execution, and output handling required for each sample.

Basic usage
============

There are typically three main steps involved in the process of generating datasets: initializing the problem class, adding various parameters and running simulation for different parameter values. To get started with the ``Airfoil`` problem, you will need airfoil coordinates stored in a ``dat`` file. There are few important points to note regarding the ``dat`` file:

- The ``dat`` file must follow the selig format i.e. the points should start from trailing edge and go in counter-clockwise direction, and then back to trailing edge.

- The first and last point in the ``dat`` file should be same (for both sharp and blunt trailing edge). This ensures that airfoil surface created using those points is closed. It is recommended to have a blunt trailing edge.

- The number of points and their distribution defines the surface mesh for the airfoil. So, if you want specific distribution of points (e.g. more points near leading edge), then make sure that the points are distributed accordingly.

- The leading and trailing edge of the airfoil coordinates must lie on the x-axis.

- All the x-coordinates must be in the range of [0,1].

All the files used in below sections can be found in ``examples`` directory on github.

Initializing the problem
-------------------------

The following code snippet shows an example of how to initialize an ``Airfoil`` problem:

.. literalinclude:: ../../../examples/airfoil_adflow_cst/runscript.py
    :start-after: # rst INIT start
    :end-before: # rst INIT end

Let's go through the code step-by-step.

First, a ``solver_options`` dictionary is created to specify the settings used by ``ADflow``. The appropriate settings depend on the specific case being simulated, so make sure to set these options appropriately. Refer to `ADflow options <https://mdolab-adflow.readthedocs-hosted.com/en/latest/options.html>`_ for more details.

Next, an ``AeroProblem`` object is created to define the baseline flow conditions for the simulation. In this example, the angle of attack, Mach number, Reynolds number, temperature, and reference quantities are specified. ``AeroProblem`` is provided by the `MDOLab Baseclasses <https://mdolab-baseclasses.readthedocs-hosted.com/en/latest/pyAero_problem.html>`_ package and supports several different ways of specifying the flow conditions. Refer to the ``AeroProblem`` `documentation <https://mdolab-baseclasses.readthedocs-hosted.com/en/latest/pyAero_problem.html>`_ for a complete description of the available options.

Next, the ``AirfoilADflowCSTOptions`` dataclass is initialized which provides an interface for defining various configuration options for the ``Airfoil`` problem. This dataclass consists of various pre-defined options, along with three mandatory options: ``airfoil_file``, ``solver_options``, and ``aero_problem``. Refer :ref:`here <problems/airfoil:Options>` for a complete description of the available options. Finally, an ``AirfoilADflowCST`` object is created using the configured options. This object represents the initialized airfoil problem that provides an interface for defining various parameters and for evaluating different parameter sets.

The general initialization workflow can therefore be summarized as:

#. Define the ADflow solver settings using solver_options.
#. Define the baseline flow conditions using AeroProblem.
#. Creat the ``AirfoilADflowCSTOptions`` dataclass object using appropriate options.
#. Create the ``AirfoilADflowCST`` problem using the configured dataclass object.

Adding parameters
------------------

Next step is to add parameters that can be varied to generate datasets. The ``add_parameter`` method is used for adding parameters to the ``Airfoil`` problem. The ``add_parameter`` method requires three arguments:

- ``name (str)``: name of the parameter to be added. The possible parameter names are: 

    - ``lower_cst``: lower surface cst coefficients, the ``num_cst_lower`` option governs the number of CST coefficients

    - ``upper_cst``: upper surface cst coefficients, the ``num_cst_upper`` option governs the number of CST coefficients

    - ``alpha``: angle of attack of the flow

    - ``mach``: mach number of the flow

    - ``reynolds``: reynolds number for the flow

- ``lower_bound (numpy array or float)``: lower bound for the variable

- ``upper_bound (numpy array or float)``: upper bound for the variable

.. note:: To add ``mach``, ``alpha``, or ``reynolds`` as a parameter, it must be defined as one of the arguments while creating ``AeroProblem`` object
 
The following code snippet shows an example of how to add parameters:

.. literalinclude:: ../../../examples/airfoil_adflow_cst/runscript.py
    :start-after: # rst PARM start
    :end-before: # rst PARM end

Evaluating samples
-------------------

After adding parameters, evaluating samples is straight forward. Following code snippet illustrates how to evaluate samples:

.. literalinclude:: ../../../examples/airfoil_adflow_cst/runscript.py
    :start-after: # rst RUN start
    :end-before: # rst RUN end

First, 10 samples are generated using a latin hypercube sampling (LHS) method within a unit hypercube. These samples are then scaled to the correct bound for the problem. The bounds for a problem can be accessed using the ``bounds`` property method from the initialized object. It returns a tuple containing two 1-D numpy arrays: first array is the lower bound and second array is the upper bound. 

In the code snippet shown above, ``samples`` is a numpy array of shape ``(10,15)``. The ``10`` indicates the number of samples while ``15`` denotes the number of parameters. Essentially, each row in the ``samples`` array is a different sample, i.e., a different set of parameter values. The order of values in each sample array depends on the order in which parameters are added. For example, in this case, the first six entries correspond to lower CST parameter, next six entries are for upper CST parameter, and the remaining three entries are for mach, alpha, and reynolds.

Next, the initialized ``Airfoil`` problem object (``airfoil`` in this case) is be called with the samples that are to be evaluated. The call method for ``Airfoil`` problem accepts two arguments:

- ``x (numpy array)``:  an array representing samples that will be evaluated

- ``return_results (bool)``: flag to determine if the results should be returned or not (default = ``False``). This should be set to ``True`` only when you want this function to return the results (e.g. in case of active learning)

When the 

If the ``return_results`` argument is set to ``False`` (default behaviour), then all the requested output generated from the simulation is stored in the respective folder.

Options
========

Following is the exhaustive list of options that can be set by the user for the ``Airfoil`` problem:

.. pydantic-options:: blackbox.airfoil_adflow_cst_problem.AirfoilADflowCSTOptions
