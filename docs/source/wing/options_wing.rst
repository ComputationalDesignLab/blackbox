.. _options_wing:

******************
Options
******************

Following is the list of options for ``WingFFD`` module. Options are divided
into two categories: mandatory and optional.

Mandatory arguments
--------------------

  ``solverOptions (dict)``
    options for the flow solver, either adflow or dafoam

  ``gridFile (str)``
    filename containing volume mesh of the wing in cgns format

  ``liftIndex (int)``
    an int describing lift direction. Only two possible values: when the
    value is 2, lift direction is y. When the value is 3, lift direction is z

  ``aeroProblem``
    aero-problem from baseclasses package defining information related aerodynamic analysis

Optional arguments
-------------------

  ``solver (str, default="adflow")``
    name of the solver, only two values are possible: "adflow" or "dafoam"

  ``openfoamDir (str, default=".")``
    name of the directory containing supporting openfoam folders such as 0, constant, and system.
    The default value is the current directory

  ``ffdFile (str, default=None)``
    filename containing coordinates of the ffd points in plot3d format

  ``directory (str, default="output")``
    name of the directory where the results will saved

  ``noOfProcessors (int, default=4)``
    desired number of processors to run the analysis on

  ``sliceLocation (list, default=[])``
    a list containing spanwise locations where the slice is required, only works with adflow

  ``writeDeformedFFD (bool, default=False)``
    option to specify whether to write deformed FFD file or not

  ``computeVolume (bool, default=False)``
    option to specify whether to compute volume of the wing or not. Setting this option to
    ``True`` without providing setting ``ffdFile`` option will raise an error. Note that the
    volume computed will be scaled by the original value.

  ``leList (list, default=None)``
    list of points (x,y,z) close to the LE of the wing but within the wing surface. Refer `this image <https://mdolab-mach-aero.readthedocs-hosted.com/en/latest/_images/opt_thickness_and_vol_diagram.png>`_
    to understand how the leList and teList options are used to compute volume. 

  ``teList (list, default=None)``
    list of points (x,y,z) close to the TE of the wing but within the wing surface. Refer `this image <https://mdolab-mach-aero.readthedocs-hosted.com/en/latest/_images/opt_thickness_and_vol_diagram.png>`_
    to understand how the leList and teList options are used to compute volume. 
    
  ``getSurfaceForces (bool, default=False)``
    option to specify whether to get force at each surface mesh coordinate or not

  ``alpha (str, default="explicit")``
    option to specify whether to consider alpha as an explicit or implicit variable. There are only two possbile values:
    ``explicit`` (normal analysis) and ``implicit`` (internal root finder). When this option is set to implicit, then for each sample secant method
    is used to find alpha such that target CL is achived. So, this option also takes longer to evaluate. **Note**: When this option is set to implicit, then 
    ``alpha`` cannot be added as a DV. This option is only for adflow

  ``targetCL (float, default=0.824)``
    this option only applies when ``alpha`` is set to ``implicit``. 
    This option specifies the value of target CL to be met. This option is only for adflow

  ``targetCLTol (float, default=1e-4)``
    this option only applies when ``alpha`` is set to ``implicit``. 
    This option specifies the required tolerance to be met for target CL. This option is only for adflow.

  ``startingAlpha (float, default=2.5)``
    this option only applies when ``alpha`` is set to ``implicit``. 
    This option specifies the initial guess for angle of attack (in degrees) 
    for secant method. Note that this value has a huge impact on the convergence speed of the secant method. 
    So, it is recommended to set this value close to the expected value. This option is only for adflow

  ``samplingCriterion (str, default="cm")``
    this option decides which method to use to generate Latin Hypercube samples. Only four options are available:
    ``c``, ``m``, ``cm``, ``ese``

  ``randomState (int, default=None)``
    this option is used to set the random state while generating Latin Hypercube samples
