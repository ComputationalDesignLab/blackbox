**************************
Wing sample generation
**************************

Blackbox also provides a module for generating wing samples. The ``WingFFD`` module provides various
methods for running analysis and generating wing data with FFD Parametrization. This module
supports two types of solver: `ADflow <https://github.com/mdolab/adflow>`_ and `DAFoam <https://github.com/mdolab/dafoam>`_. 
To use this module, two files are requried:

- Volume mesh in cgns format which defines computational mesh for the analysis
- FFD file in plot3d format which defines a FFD box around the wing

Unlike airfoil sample generation, whenever shape of the wing is modified using FFD, the volume mesh is deformed using
`idwarp <https://github.com/mdolab/idwarp>`_. So, ensure that idwarp is installed before using this module,
refer :ref:`installation<installation>` details.

There are other required files depending on the solver, refer below sections to know more about 
how to use this module. All the examples used in below sections can be found in ``examples`` directory on github.

.. toctree::
    :maxdepth: 3
    :caption: Table of Contents
    
    adflow
    dafoam
    options_wing