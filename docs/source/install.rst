.. _install:

Installation
============

Follow below steps for installing Blackbox:

- Clone or download the latest tagged release from Blackbox's `github repository <https://github.com/ComputationalDesignLab/blackbox/releases>`_.
- Open the terminal and ``cd`` into the root of cloned/downloaded repository
- Activate the virtual environment in which you want to install Blackbox
- Run the following command to instal the package::

	pip install .

- If you want to install the package in editable mode, run the following command::
	
	pip install -e .

.. note:: This will not install additional packages/libraries required for specific problems, you will need to install those separately. Refer specific problem section for more details about required packages/libraries.
