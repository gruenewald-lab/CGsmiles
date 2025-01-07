Layering of Resolutions
=======================

CGSmiles enables the representation of molecular graphs at arbitrary resolutions
and their connection to progressively finer resolutions, allowing for the
hierarchical layering of multiple levels of details.

Basic Syntax Features
---------------------

Base Graph & Resolutions
^^^^^^^^^^^^^^^^^^^^^^^^
The notation starts with the coarsest representation of the system – the base
graph. This graph is enclosed in curly braces. Each additional resolution is
represented as a list of fragment graphs, also enclosed in curly braces and
separated from the preceding resolution graph by a period. If the final resolution
graph is at the atomic level, either CGSmiles or OpenSMILES syntax can be used
to describe the fragment graph. This dual approach allows seamless conversion to
atomistic resolution using established standards, while also supporting
intermediate coarse-grained representations.

.. code-block:: none

  {coarsest-graph}.{fragments-resolution-1}.{fragments-resolution-2}

Advanced Syntax Features
------------------------

Linearizing Rings
^^^^^^^^^^^^^^^^^
Rings at the atomistic resolution can often be mapped into linear structures
at the CG level, a common practice in chemically specific force fields such
as Martini. In the CGSmiles notation, bond orders at the coarser resolution are
utilized to describe such a case.

For example, cyclohexane is represented at the Martini 3 level with a bond
order of 2. This indicates that at the next finer resolution level, two bonds
must connect the atoms corresponding to the two CG nodes.

.. code-block:: none

  Martini3 Cyclohexane
  {[#SC3]=[#SC3]}.{#SC3=[$]CCC[$]}

This approach also extends to more complex cases, such as  splitting fused rings
with three or more shared bonds at the CG level. Each additional ring increases
the bond order.

.. code-block:: none

  Napthalene split as two particles
  {[#A]#[#A]}.{#A=[$]CCC[$]CC[$]}

The current CGSimles syntax supports bond orders up to 4, which defines the
maximum number of ring connections that can be represented linearly.

Virtual Edges
^^^^^^^^^^^^^
In certain scenarios, a CG model might include interacting particles that do not
correspond to any finer-resolution nodes or atoms. For example, at the Martini 3
resolution glucose is represented by three CG particles splitting the sugar ring
and one additional virtual particle. The TC4 bead captures the hydrophobic
interactions at the ring center but lacks any corresponding fragments at finer
resolution. To accommodate such particles, the CGSmiles notation employs zero
bond order edges, referred to as virtual edges.

.. code-block:: none

  Martini 3 Glucose
  {[#SP4r]1.2[#SP4r].3[#SP1r]1.[#TC4]23}.{#SP4r=OC[$]C[$]O,#SP1r=[$]OC[$]CO}

Virtual edges are ignored when establishing connections and any particle with only
virtual edges is excluded entirely when transitioning to finer resolutions. We
note that these virtual edges and virtual particles are not to be confused with
the GROMACS virtual sites. A virtual site in GROMACS describes how a particle's
coordinates are constructed. If a virtual side describes real atoms or CG
particles they would be treated as regular nodes rather than virtual ones.

Overloading Wildcards
^^^^^^^^^^^^^^^^^^^^^
In certain cases, a single CG graph might describe more than one molecule at
the fine-grained resolution because of a loss in resolution at the CG level.
An example are Martini lipids such as POPC. POPC can describe lipids with a
tail length of 16 or 18 carbons and thus represents at least four molecules
when accounting for the position for the double bond. To capture this feature
CGSmiles allows to overload the wildcard (*) syntax using annotations. In
OpenSMILES a wildcard means any atom can be placed at the wildcard position.
To specify a selection of atoms CGSmiles allows to annotate a wildcard using the
select keyword abbreviated as ‘s’. Thus, a tail bead in POPC could be written as
``C1=CCCC[*;s=C,0][*;s=C,0]``. Note that the current molecule resolver is not able
to handle wildcard overloading.
