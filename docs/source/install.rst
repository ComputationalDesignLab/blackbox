.. _install:

Installation
============

Blackbox uses various packages for sample generation and many of them are only supported on Linux-based systems.
Due to this dependency, Blackbox is limited to these systems only. The dependencies are divided into two categories:

- **MACH-Aero packages**: `MACH-Aero <https://mdolab-mach-aero.readthedocs-hosted.com/en/latest/index.html>`_ is a 
  framework for performing aerodynamic shape optimization. It consists of various packages, some of which are used
  by Blackbox while generating samples. Below table outlines the requried packages:

  .. table::
    :widths: auto

    =========================================================== ================
    Package                                                     Version
    =========================================================== ================
    `baseclasses <https://github.com/mdolab/baseclasses>`_      1.7.0 or higher
    `pyGeo <https://github.com/mdolab/pygeo>`_                  1.12.2 or higher
    `pyHyp <https://github.com/mdolab/pyhyp>`_                  2.6.0 or higher
    `idwarp <https://github.com/mdolab/idwarp>`_                2.6.1 or higher
    `cgnsutilities <https://github.com/mdolab/cgnsutilities>`_  2.7.1 or higher
    `ADflow <https://github.com/mdolab/adflow>`_                2.7.4 or higher
    =========================================================== ================

  Please refer MACH-Aero's `installation page <https://mdolab-mach-aero.readthedocs-hosted.com/en/latest/installInstructions/installFromScratch.html>`_ 
  for more details about installing these packages.

- **Other packages**: Other than MACH-Aero packages, Blackbox also requries following packages:

  .. table::
    :widths: auto

    ==================================================== ================
    Package                                              Version
    ==================================================== ================
    `prefoil <https://github.com/mdolab/prefoil>`_       2.0.0
    `smt <https://smt.readthedocs.io/en/stable/>`_       2.0.0 or higher
    `psutil <https://psutil.readthedocs.io/en/latest/>`_ 5.9.1 or higher
    ==================================================== ================

It is recommended to create a virtual environment for installing all the above packages. Follow below steps for installing Blackbox:

- Clone or download the latest tagged release from Blackbox's `github repository <https://github.com/ComputationalDesignLab/blackbox>`_.
- Open the terminal and ``cd`` into the root of cloned/downloaded repository
- Activate the virtual environment and run::

    pip install .
