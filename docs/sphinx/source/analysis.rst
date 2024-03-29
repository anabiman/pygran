********
Analysis
********

Submodules
##########
.. autosummary::

   :toctree: _autosummary

   pygran_analysis.core
   pygran_analysis.dynamics
   pygran_analysis.equilibrium


Introduction
############

The `analysis <https://github.com/Andrew-AbiMansour/PyGranAnalysis>`_ module enables programmers to read and analyze DEM trajectory files. The trajectory file stores sequential snapshots of the system positions, velocities, forces, stresses, or other physical quantities of interest. The most fundamental class in this module is `System <autosummary/analysis.core.html#analysis.core.System>`_, which uses a factory class to instantiate a `SubSystem <autosummary/analysis.PyGran.core.html#analysis.core.SubSystem>`_ subclass. In principal, any implementation of `SubSystem`_ can be instantiated with this factory. A `System`_ instance is always created by passing the filename of a specific `SubSystem`_::

	Granular = System(SubSystem='path/to/file')

`System`_ is an iterator. Thus, it can be iterated/looped over when the supplied `SubSystem`_ contains a time series. For example, the statement::

	for frame in Granular:
		# do something with Granular.SubSystem.property

loops over every frame stored in the suppled trajectory file, returning the frame number at each instant.

`SubSystem`_ is an abstract class that encapsulates common attributes and methods for basic DEM objects such as `Particles <autosummary/analysis.core.html#analysis.core.Particles>`_ and `Mesh <autosummary/analysis.core.html#analysis.core.Mesh>`_, both being derived classes of `SubSystem`_. While `SubSystem`_ is a mutable object, its properties cannot be directly modified by the user, i.e. they can modified only by the methods in `SubSystem`_. The basic data structure in this object is a Python dictionary (``SubSystem.data``) which contains references to the `SubSystem`_  attributes and is used to generate the dynamic interface of a `SubSystem`_ object. The attributes of `SubSystem`_  change from one frame to another for the same system, thus, ``SubSystem.data`` is updated every single time the `System`_ is evolved in time.

Systems and SubSystems form the core classes used in the `analysis`_ module (:numref:`fig_classes`), and they provide numerous methods for temporal, structural, and image analyses as covered in more detail in the next sections.

.. _fig_classes:
.. figure:: images/classes.png
    :scale: 64%
    :align: center
    :alt: basic classes of PyGran.analysis module
    :figclass: align-center

    A UML diagram of the three fundamental classes and some of their methods and attributes in the analysis module. The state of a DEM system is determined by the System class, which creates and stores one or more instances of classes derived from SubSystem that describe basic DEM elements such as particles or surface mesh(es).

Systems & SubSystems
####################

System constructor
~~~~~~~~~~~~~~~~~~
The `System`_ class is the most fundamental class in *PyGran*. It uses a factory to create objects derived from `SubSystem`_  that describe the state of a granular system (:numref:`fig_classes`). These subclasses can be instantiated from an input *data* dictionary or copied from another instance of `SubSystem`_.
`System`_ creates an instance (or a list of instances) of `SubSystem`_ from input filename strings (or list of strings) that are passed to ``System.__init__`` by a factory object.

 `System`_ contains all the objects, methods, and properties that describe the state of a DEM system. This class also handles I/O operations and ensures proper frame to frame  propagation when reading input trajectory files. The frame is controlled only by `System`_ when the latter is looped over via methods defined in a `SubSystem`_ sublass (read/write functions). Since DEM simulations consist of a set of particles in contact with surface triangulations (representing walls), `System`_  creates subclasses of `SubSystem`_ such as `Particles`_  and \emph{Mesh} (:numref:`fig_classes`) based on input trajectory files. The 4 different unit systems supported by this class are summarized in Table :numref:`table_units`.

SubSystem constructor
~~~~~~~~~~~~~~~~~~~~~
This is an abstract class that encapsulates common attributes and methods for basic DEM objects such as `Particles`_  and `Mesh`_, both being derived classes of `SubSystem`_. While `SubSystem`_ is a mutable object, its properties cannot be directly modified by the user, i.e. they can modified only by the methods in `SubSystem`_. The basic data structure in this object is a Python dictionary (``SubSystem.data``) which contains the `SubSystem`_ attributes and is used to instantiate a `SubSystem`_ object, i.e. ::

	NewSS = SubSystem(**input_data)

Alternatively, `SubSystem`_ objects can be used to create new `SubSystem`_ objects (i.e. copy constructor)::

	CopySS = SubSystem(SubSystem=OriginalSS)

.. todo::
	A `System`_ object can also be sliced (by frames), e.g. the following statement ::
	
		SlicedSys = System[start:end]

	yields a new `System`_ object (SlicedSys) that contains all frames from *start* to *end-1*.

Particles
~~~~~~~~~
The `Particles`_ class provides a way to store, manipulate, and operate on particle attributes generated by DEM simulation. This class is a subclass of  
`System`_ and can therefore be sliced and looped over. Furthermore, this class provides several basic routines for computing properties usually encountered in powder technology (such as mass density, radial distribution function, radius of gyration, etc.) as well as particle-based operators discussed below.

Binary operations
~~~~~~~~~~~~~~~~~
Extended assignments can be made to `Particles`_  with ``+=``. For example, *Particles_i* is appended to ``Particles``  with the following statement::

	Particles += Particles_i

If ``Particles_i`` has fewer attributes than those in ``Particles`` , then this assignment is rejected. Otherwise, any additional attributes of ``Particles_i`` not found in ``Particles``  are neglected.

2 `Particles`_  objects can be concatenated with the ``+`` operator. This operation can lead to reduction in the number of attributes if one of the classes being added has fewer attributes than the other(s). In this case, the resultant `Particles`_  will acquire concentenated attributes specified by the class with minimum number of attributes.  2 `Particles`_  objects can also be multiplied wth ``*`` to yield a new object whose vector attributes are the geometric mean of the external product of the vector attributes of the two objects being multiplied. For instance, if 3 objects ``Particles_i``, ``Particles_j``, and ``Particles_k`` contain :math:`n_i`, :math:`n_j`, and :math:`n_k` particles, respectively, then the following code ::

	Particles = Particles_i + Particles_j * Particles_k

yields a new `Particles`_  object containing :math:`n_i + n_j n_k` particles and with vector attributes :math:`[a_{i,1}, ... , a_{i,n_i}, \sqrt{a_{j,1} \times a_{k,1}}, ... \sqrt{a_{j,n_j n_k} \times a_{k,n_j n_k}}]`.

Basic methods
~~~~~~~~~~~~~
Some of the basic methods available to `Particles`_  are shown in :numref:`fig_classes`. Furthermore, the ``PyGran.analysis`` module provides a `Neighbors <autosummary/analysis.core.html#analysis.equilibrium.Neighbors>`_ class that is instantiated with a `Particles`_  object to provide methods for nearest neighbor analysis. With this class, properties such as coordination numbers, overlap distances, and force chains can be readily computed.

Input/output
************
Any class derived from `SubSystem`_ must implement read/write methods. In the current version, *PyGran* supports reading and writing particle trajectory files for *LIGGGHTS*. The input trajectory can be a dump or a vtk :cite:`schroeder2004visualization` file.

Custom SubSystems
*****************
User-defined subclasses of `SubSystem`_ can be easily created by using Python's inheritence feature. The keyword ``module`` must be passed to the subclass constructor in order to make sure *PyGran* imports the module containing the subclass.

*PyGran*'s extensible and object-oriented design makes it ideal for creating user-defined particles. Since `System`_ uses a Factory class to instantiate a `Particles`_  or `Mesh`_ object, it can in principle be used to instantiate a user-defined class. This is demonstrated in the code below for a simple coarse-grained class that demonstrates the use of the ``filter`` method to eliminate particles overlapping by a certain %.

A simple user-defined ``CoarseParticles`` class can be defined as a subclass of `Particles`_  with two key arguments: ``scale``, which controls the level of coarse-graining (or reduction) and ``percent`` which is used to eliminate the resultant coarse-grained particles overlapping by a certain percentage with respect to their radius. A script that implements this class is shown below::

	class CoarseParticles(analysis.Particles):
        	def __init__(self, **args):
                	super().__init__(**args)

                	if 'scale' in args and 'percent' in args:
                        	self.scale(args['scale'], ('radius',))
                        	CG = analysis.equilibrium.Neighbors(self).filter(percent=args['percent'])

                        	self.__init__(CoarseParticles=CG)

The ``CoarseParticles`` object uses a recursive call to instantiate a derivative of the `Particles`_  class and therefore inherits all of the latter's properties and methods.

Mesh
====
The `Mesh`_ class uses the VTK library to read input mesh files and expose the stored attributes (nodes, positions, stresses, etc.) to the user.

Surface walls are represented in *PyGran* by the `Mesh`_ class, a subclass of `SubSystem`_ (see :numref:`fig_classes`). This class uses the VTK library :cite:`schroeder2004visualization` to read an input mesh trajectory (one or more sequence of VTK file(s)) and expose all of the stored file variables to the user. This is particularly useful for analyzing DEM simulation involving mesh-particle interaction or coupled CFD-DEM simulations.

