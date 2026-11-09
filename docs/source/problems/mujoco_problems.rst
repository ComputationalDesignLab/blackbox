================
MuJoCo problems
================

MuJoCo (Multi-Joint dynamics with Contact) is a physics engine for simulating robotics and other systems with articulated bodies and contact interactions. Blackbox uses MuJoCo through `Gymnasium <https://gymnasium.farama.org/environments/mujoco/>`_, which provides a consistent API for creating an environment, resetting it, and stepping through a simulation. The following Gymnasium MuJoCo environments are supported within Blackbox:

* `Half-Cheetah <https://gymnasium.farama.org/environments/mujoco/half_cheetah/>`_
  (``HalfCheetah-v5``)
* `Hopper <https://gymnasium.farama.org/environments/mujoco/hopper/>`_
  (``Hopper-v5``)
* `Humanoid <https://gymnasium.farama.org/environments/mujoco/humanoid/>`_
  (``Humanoid-v5``)
* `Swimmer <https://gymnasium.farama.org/environments/mujoco/swimmer/>`_
  (``Swimmer-v5``)
* `Walker2d <https://gymnasium.farama.org/environments/mujoco/walker2d/>`_
  (``Walker2d-v5``)

You must install Gymnasium and its MuJoCo dependencies separately by running following command in terminal::

  pip install gymnasium[mujoco]

Basic usage
============

The example below demonstrates how to work with ``MujocoProblem`` class implemented in Blackbox.

.. literalinclude:: ../../../examples/mujoco_problems/half_cheetah.py
    :start-after: # rst start
    :end-before: # rst end

The problem class is initialized with the environment name, number of rollouts (number of times the experiment is repeated), and a seed. The ``env_name`` must be one of the following: ``HalfCheetah-v5``, ``Hopper-v5``, ``Humanoid-v5``, ``Swimmer-v5``, ``Walker2d-v5``. There are other arguments as well that can be provided during initialization, refer to mujoco problem :class:`API <blackbox.mujoco_problems.MujocoProblem>` for more details. Then, a Latin hypercube sampling plan is created and scaled to appropriate bounds. In this example, the ``x`` is of shape (10,102), where 10 corresponds to number of samples and 102 denotes problem dimension. Lastly, the sampling plan is evaluated by calling the initialized class object. The ``call`` method can takes a numpy array of shape ``(dim,)`` (single sample) or ``(n_samples, dim)`` (multiple samples), and returns the value of reward function associated with each sample. You can save the ``x`` and ``y`` as an ``npz`` or ``mat`` file for later use. Similar to this example, you can initialize mujoco problem with any of the valid environments and generate datasets.