===================
Aerostruct Problem
===================

The ``Aerostruct`` problem provides an interface for generating static aerostructural wing datasets using a coupled high-fidelity analysis. Specifically, this problem couples a computational fluid dynamics (CFD) solver with a finite element (FE) solver to perform static aerostructural analysis. The ``Aerostruct`` problem currently uses ADflow_ as the CFD solver, TACS_ as the FE solver, and OpenVSP for wing shape parameteriation. Users can parameterize wing planform parameters, wing airfoil cross-section, operating conditions, and evaluate the resulting configurations without manually setting up the geometry paramerization, coupled analysis, and post-processing workflow. Additional CFD solvers, such as OpenFOAM and SU2, may be supported in future releases. The following sections describe the dependencies, computational workflow, and how to generate aerostructural datasets with the ``Aerostruct`` problem.

Dependencies
=============

The ``Aerostruct`` problem requires several additional packages that must be installed before using this problem. The required packages and their corresponding versions are listed below.

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Package
     - Version
     - Comments
   * - `mpi4py`_
     - 3.1.6
     - Used for parallel processing during mesh deformation, CFD, and FE solver execution
   * - `psutil`_
     - `-`
     - Used for monitoring child processes launched for coupled analysis
   * - `mdolab-baseclasses`_
     - >= 1.8.4
     - Provides ``AeroProblem`` for defining flow conditions
   * - `openvsp`_
     - `-`
     - Handles wing shape paramerization
   * - `pygeo`_
     - >= 1.17.0
     - Provides an interface for handling wing shape parameterization through OpenVSP
   * - `idwarp`_
     - >= 2.6.2
     - Used for morphing volume mesh based on given shape parameters and during coupled analysis
   * - `adflow`_
     - 2.11.0
     - Used as the CFD solver for aerodynamic simulations
   * - `tacs`_
     - 3.7.2
     - Used as the FE solver for structural simulations
   * - `funtofem`_
     - 0.3.7
     - Used for transfering loads and displacement between CFD/FE solvers
   * - `mphys`_
     - 1.3.0
     - Used for setting up coupled aerostructural analysis with CFD and FE solver
   * - `pyvista`_
     - `-`
     - Used for reading, processing, and writing surface and volume field data
   * - `h5py`_
     - `-`
     - Used for creating and storing field results in HDF5 format

.. _mpi4py: https://mpi4py.readthedocs.io/en/3.1.6/
.. _mdolab-baseclasses: https://github.com/mdolab/baseclasses
.. _openvsp: https://github.com/OpenVSP/OpenVSP/tree/main/src/python_api/packages#openvsp-python-packages
.. _pygeo: https://github.com/mdolab/pygeo
.. _idwarp: https://github.com/mdolab/idwarp
.. _ADflow: https://github.com/mdolab/adflow
.. _TACS: https://github.com/smdogroup/tacs
.. _pyvista: https://docs.pyvista.org/
.. _h5py: https://docs.h5py.org/
.. _funtofem: https://github.com/smdogroup/funtofem
.. _mphys: https://github.com/openMDAO/MPHYS
.. _psutil: https://pypi.org/project/psutil/

Make sure that all required packages are installed with compatible versions before using the ``Aerostruct`` problem.

Workflow
=========

The ``Aerostruct`` problem accepts a set of samples, where each sample contains parameter values. This set may contain a single sample or multiple samples. For each sample :math:`x_i` in the sample set, the ``Aerostruct`` problem performs the sequence of operations described below:

1. **Wing shape parameterization**: Geometry parameters are used to define the wing planform shape and airfoil cross-section through ``OpenVSP`` and ``pyGeo``. Based on the geometry parameters, baseline CFD surface mesh coordinates and baseline wingbox FE mesh are deformed.
2. **Volume Mesh deformation**: The baseline CFD volume mesh is deformed based on the updated surface mesh coordinates using ``idwarp``.
3. **Coupled analysis**: A coupled analysis is then performed using the updated CFD volume mesh and the wingbox FE mesh, along with flow parameters and load factor.
4. **Output extraction**: Scalar and field quantities from the converged solver results are extracted and stored in a folder.

When multiple samples are provided, the workflow is evaluated for each sample, and the resulting outputs are stored separately for each sample. This allows the ``Aerostruct`` problem to provide a consistent interface for generating datasets without requiring users to manually manage geometry processing, mesh deformation, coupled solver execution, or output handling.

The complete workflow is illustrated in the XDSM diagram below.

.. image:: ../_static/xdsm_aerostruct.png
    :width: 100%
    :align: center
    :alt: aerstruct problem xdsm
    :class: image-margin

The stacked rectangles representing various blocks such as shape paramerization, volume mesh deformation, CFD, and FE analysis indicate that these operations use multiple processors through MPI. The ``Aerostruct`` problem therefore abstracts geometry processing, mesh deformation, coupled-solver execution, and output handling required for each sample.

Basic usage
============

Setting up ``Aerostruct`` problem and generating dataset is very similar to ``Wing`` problem. It is recommended that you go through ``Wing`` problem :ref:`basic usage <wing-basic-usage>` section before moving forward. To use the ``Aerostruct`` problem, you will require following items:

1. Baseline volume CFD mesh in CGNS format
2. Baseline wingbox FE mesh in BDF format
3. OpenVSP model for the baseline wing (only required if wing shape will be parametrized)

Make sure that the CFD and FE mesh work with corresponding solvers before using it with ``Aerostruct`` problem. The requirements for the OpenVSP file are provided :ref:`here <wing-openvsp-guidelines>`. The files used in the example demonstrated below are available in the ``examples/wing_aerostruct_vsp`` directory.

Defining flow and wing geometry parameters
-------------------------------------------

The flow and wing geometry parameters are defined in a similar manner as it is done for ``Wing`` problem:

* :ref:`Flow parameters <wing-flow-parameters>`: create an ``AeroProblem`` object and add flow parameters
* :ref:`Wing geometry parameters <wing-geometry-parameters>` : create ``WingVSP`` object and add wing geometry parameters

.. note:: Load factor can be added as a parameter while initializing the ``Aerostruct`` problem
 
Following code snippet illustrates how to define flow and wing geometry parameters:

.. literalinclude:: ../../../examples/wing_aerostruct_vsp/runscript.py
    :start-after: # rst PARAM start
    :end-before: # rst PARAM end

Initializing the problem
-------------------------

Next step is to initialize the ``WingAeroStructOptions`` class by setting appropraite options. There are several options that can be provided while creating the options class, refer to :ref:`complete options list <problems/aerostruct:Options>` for all available options. Most of the options are optional, except the ``aero_mesh_file``, ``aero_solver_options``, ``aero_problem``, and ``struct_mesh_file`` are mandatory options. The initialized options object is used to create the main ``WingADFflowVSP`` problem object. Following code snippet illustrates how to initialize the ``WingAeroStructVSP`` problem class:

.. literalinclude:: ../../../examples/wing_aerostruct_vsp/runscript.py
    :start-after: # rst INIT start
    :end-before: # rst INIT end

Note that the load factor is added as a parameter by setting ``load_factor_parameter`` to ``True``. You can also set its lower and upper bound using ``load_factor_lb`` and ``load_factor_ub``, respectively.

Evaluating samples
-------------------

After defining the problem, evaluating samples is straightforward. Following code snippet illustrates how to evaluate a set of samples:

.. literalinclude:: ../../../examples/wing_aerostruct_vsp/runscript.py
    :start-after: # rst RUN start
    :end-before: # rst RUN end

Firstly, 10 samples are generated using a Latin hypercube sampling (LHS) method within a unit hypercube. You can use any sampling method of your choice to generate samples. These samples are then scaled to appropriate bounds for the problem. The bounds can be accessed through the ``bounds`` property of the initialized problem object. This property is a tuple containing two 1D numpy arrays: first array contains the lower bounds and the second array contains the upper bounds.

In the code snippet above, samples is a numpy array with shape (10, 20). The 10 represents the number of samples, while 20 represents the number of parameters. Each row in samples array corresponds to a different sample, i.e., a different set of parameter values. The order of the values in each row is determined by the order in which the parameters were added to the problem.

For example, in this case, the first three entries correspond to flow parameters. The next four entries are for wing planform parameters. The next twelve entries correspond to cross-section airfoil CST parameters, and the remaining entry corresponds to load factor. Note that, in aerostruct problem, flow parameters will always preceed shape parameters, while the load factor will always be the last parameter, if added as a parameter.

The initialized ``Aerostruct`` problem object (wing_aerostruct in this case) can then be called with the samples to be evaluated. The ``__call__`` method accepts two arguments:

- ``x (numpy array)``: an array of shape (N,D) or (D,) representing the samples to be evaluated. Here, ``N`` denotes the number of samples and ``D`` denotes the number of parameters

- ``return_results (bool)``: a flag indicating whether the results should be returned by the function (default = ``False``). Set this to ``True`` when the results are needed directly in Python, such as during active learning. When ``return_results`` is ``False``, the requested outputs from each simulation are stored in the corresponding sample directory.

Once the dataset generation process starts, a directory with the name specified by the ``directory`` option is created. All generated data is stored within this directory. The directory structure is organized as follows:

- The main directory contains a ``description.txt`` file that provides information about the generated dataset, including the parameters and their bounds.

- A separate subdirectory is created for each evaluated sample. The subdirectories are numbered sequentially starting from ``1`` (i.e., ``1``, ``2``, ``3``, and so on), with each number corresponding to the row of the samples array that was evaluated. For example, the directory ``1`` contains the results for ``samples[0]``, ``2`` contains the results for ``samples[1]``, and so on.

- Each sample directory contains the files generated during the simulation. Depending on the options provided during initializating ``WingAeroStructOptions`` class, the following files may be present:

  - ``log.txt``: a log file containing messages generated during mesh deformation and solver execution
  - ``parameters.json``: a JSON file containing the parameter values used for the sample
  - ``scalar_outputs.json``: a JSON file containing the scalar outputs specified by the ``scalar_outputs`` option.
  - ``surface_solution.cgns``: the surface solution written by the CFD solver. This file is stored only when ``write_surface_output`` is set to ``True``
  - ``volume_solution.cgns``: the volume solution written by the CFD solver. This file is stored only when ``write_volume_output`` is set to ``True``
  - ``structural_solution.f5``: the structural solution written by the FE solver. This file is stored only when ``write_structural_output`` is set to ``True``
  - ``surface_outputs.hdf5``: field data extracted from the surface solution and stored in an HDF5 file. This file is stored only when ``write_surface_output`` is set to ``True``. It will contain entities provided in the ``surface_outputs`` option
  - ``volume_outputs.hdf5``: field data extracted from the volume solution and stored in an HDF5 file. This file is stored only when ``write_volume_output`` is set to ``True``. It will contain entities provided in the ``volume_outputs`` option
  - ``structural_solution.vtu``: structural data extracted from the ``f5`` file and stored in ``vtu`` format. This file is stored only when ``write_structural_output`` is set to ``True``.
  - ``lift_distribution.dat``: a dat file containing lift distribution along the wing span, this file is written only when ``write_lift_distribution`` is set to ``True``
  - ``slices.dat``: a dat file containing flow quantities along wing cross-section at different locations (from root till tip), this file is stored only when ``write_slice_file`` is set to ``True``
  - ``updated_model.vsp3``: an OpenVSP file containing updated wing model, this file is stored only when ``write_vsp_file`` is set to ``True``
  - ``updated_model.stl``: a stl file describing the updated wing model, this file is stored only when ``write_stl_file`` is set to ``True``

Options
========

Following is the exhaustive list of options that can be set by the user for the ``AeroStruct`` problem:

.. pydantic-options:: blackbox.aerostruct_problem.WingAeroStructOptions
