**************************
Airfoil sample generation
**************************

Blackbox provides two modules for generating airfoil samples depending on the parametrization. 
The ``AirfoilCST`` and ``AirfoilFFD`` module provides various methods for running analysis and 
generating data with CST and FFD Parametrization, respectively. Before starting to use these 
modules, a .dat file is needed which contains airfoil coordinates. There are few important 
points to note regarding .dat file:

- The .dat file should follow **selig** format i.e. the points should start from trailing edge and
  go in counter-clockwise direction and then back to trailing edge.
- The first and last point in the .dat file should be same (for both sharp and blunt trailing edge). 
  This ensures that surface created using those points is closed. The coordinates obtained from
  the UIUC database usually doesn't close the loop, so check the .dat file properly before
  using with Blackbox.
- The leading and trailing edge of the airfoil coordinates should lie on the x-axis (best case scenario)
  or should be very close to it.

  .. note::
    It is highly recommended that the first and last point should be (0,0) in the .dat file. If the
    airfoil has blunt trailing edge, then the number of points along the trailing edge should be odd.

Explore below sections to know more about how these modules can be used. It is highly recommended to go through
atleast one of the first two sections before moving on to other sections. All the files used in below sections
can be found in ``examples`` directory on github.

.. toctree::
    :maxdepth: 3
    :caption: Table of Contents
    
    cst
    ffd
    single
    field_data
    imp_alpha
    options
