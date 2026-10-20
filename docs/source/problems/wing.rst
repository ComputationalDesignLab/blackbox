=================
Wing Problem
=================

The ``Wing`` problem provides an interface for generating aerodynamic datasets for wings using high-fidelity computational fluid dynamics (CFD) workflow. Users can parameterize wing planform parameters, wing airfoil cross-section, flow conditions, and evaluate the resulting configurations without manually setting up the geometry, meshing, solver, and post-processing workflow. The ``Wing`` problem currently uses ADflow_ as the CFD solver and the OpenVSP for wing shape parameteriation. Additional CFD solvers, such as OpenFOAM and SU2, may be supported in future releases. The following sections describe the dependencies, computational workflow, and how to generate aerodynamic datasets with the ``Wing`` problem.

Dependencies
=============

The ``Wing`` problem requires several additional packages that must be installed before using this problem. The required packages and their corresponding versions are listed below.

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
   * - `openvsp`_
     - `-`
     - Handles wing shape paramerization
   * - `pygeo`_
     - >= 1.17.0
     - Provides an interface for handling wing shape parameterization through OpenVSP
   * - `idwarp`_
     - >= 2.6.2
     - Used for morphing volume mesh based on given shape parameters
   * - `adflow`_
     - 2.11.0
     - Used as the CFD solver for airfoil flow simulations
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
.. _pyvista: https://docs.pyvista.org/
.. _h5py: https://docs.h5py.org/

Make sure that all required packages are installed with compatible versions before using the ``Wing`` problem.

Workflow
=========

The ``Wing`` problem accepts a set of samples, where each sample contains parameter values. This set may contain a single sample or multiple samples. For each sample :math:`x_i` in the sample set, the ``Wing`` problem performs the sequence of operations described below:

Options
========

Following is the exhaustive list of options that can be set by the user for the ``Wing`` problem:

.. pydantic-options:: blackbox.wing_problem.WingADflowVSPOptions