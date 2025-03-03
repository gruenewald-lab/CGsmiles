General Fragment Syntax
=======================

Overview
--------
CGsmiles supports the representation of molecular structures at different
resolutions through a fragment replacement syntax. This allows users to specify
more detailed molecular structures connected to a coarse graph representation.

Fragment Graphs
---------------
The notation for a fragment graph starts with a ‘#’ followed by the label of
the coarser-resolution node and an ‘=’ sign. Each fragment name must be unique
to ensure unambiguous identification. Fragments graphs use the same general
graph syntax as outlined before, however, it is permitted to use OpenSMILES
syntax to define an atomic resolution fragment.

.. code-block:: none

   Example: Bezene as single fragment
   "{#BENZ}.{#BENZ=c1ccccc1}"

Bond Operators
^^^^^^^^^^^^^^
To define how two consecutive fragments at a finer resolution are connected,
CGsmiles  builds upon the bonding connector syntax established in BigSMILES to
avoid ambiguity. Any node or atom that connects to a neighboring fragment is
followed by one of four bonding connectors (‘$’, ‘>’, ‘<’, ‘!’) enclosed in
square brackets. In addition, any operator may be combined with an alphanumeric
label to distinguish non-equivalent operators of the same type.

- **Undirected Bonding Operator $.**
  The undirected bonding operator ‘$’ connects to any other ‘$’ operator in
  connected fragments, as specified in the coarser resolution graph. An
  undirected bonding operator may be followed by an alphanumeric label,
  ensuring that only operators with matching labels are connected.

  .. code-block:: none

    Example: PEO can connect to any other PEO on the first or second carbon
    {[#PEO=[$]COC[$]}

- **Directed Bonding Operators > and <**
  A directed bonding operator can only pair with its complementary counterpart to
  ensure the correct connectivity in asymmetric fragments. These bonding operators
  can also be annotated with an alphanumeric label for further specificity.

  .. code-block:: none

    Example: Polystyrene, where the CH2 group always connects to the CH group
    {[#PS=[>]CC[<]c1ccccc1]}

- **Shared Bonding Operator !**
  To address a common scenario in CG force fields where an atom is distributed
  between two finer resolution nodes, CGsmiles introduces the shared bonding
  operator ‘!’. In the case of toluene represented at the Martini 3 level, some
  of the ring atoms are shared between the two CG beads. When two fragments are
  connected using the shared bonding operator, the atoms at the connection point
  are merged into a single atom, retaining the bonds from both fragments.

  .. code-block:: none

    Example: Martini3 Toluene where some of the carbon atoms in Toluene are
    shared between beads
    {[#SC4]1[#TC5][#TC5]1}.{#SC4=Cc(c[!])c[!],#TC5=[!]ccc[!]}

Valency
-------
CGSimles does not enforce valency rules for atoms or nodes, allowing any atom
to be followed by multiple bonding operators. In the case of all-atom fragments,
the hydrogen count is determined only after the molecule's full connection is
established. In cases where a bond of higher order needs to be represented, a
bond order symbol should be placed between the node and bonding operator in
both fragments. For example, splitting 2-pentene into two fragments results in
{[#A][#B]}.{#A=CC=[$],#B=[$]=CCC}, where the bond order symbol ‘=’ indicates a
double bond between ethane and propane fragment.


Updated Bonding Descriptors Lookup Table
----------------------------------------
This table now includes the squash descriptor, summarizing all the bonding descriptors used in CGsmiles:

+----------------+---------------------------+--------------------------------------------------------------------+
| Descriptor     | Symbol                    | Description                                                        |
+================+===========================+====================================================================+
| Indiscriminate | `[$]`                     | Connects to any matching `$` descriptor.                           |
+----------------+---------------------------+--------------------------------------------------------------------+
| Forward bond   | `[>]`                     | Must connect with a bonding descriptor of type `<`.                |
+----------------+---------------------------+--------------------------------------------------------------------+
| Backward bond  | `[<]`                     | Designed to connect with a descriptor of type `>`.                 |
+----------------+---------------------------+--------------------------------------------------------------------+
| Alphanumeric   | `[descriptor]alphanumeric`| Adds specificity to descriptors, requiring exact matches.          |
+----------------+---------------------------+--------------------------------------------------------------------+
| Squash         | `[!]`                     | Indicates overlapping mappings; connected atoms are identical.     |
+----------------+---------------------------+--------------------------------------------------------------------+
