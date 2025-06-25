Martini Mappings
================

CGsmiles can be used to define Martini mappings. Currently there a no canonical
rules on how the CGsmiles string needs to be formatted. However, we recommend to
follow the rules listed below to store all information relevant for forward
mapping, bead assignment, and backwards mapping.

Design Guidelines
-----------------

- use the bead type as fragment name in the CGsmiles string
- if there are multiple fragments that have the same bead type append letters
  A-Z (e.g. TC5 -> TC5A)
- if an atom is part of more than one fragment use the [!] bonding operator
- annotate chiral atoms in the fragment part 
- annotate cis/trans isomers in the fragment part
- annotate weights if applicable
- annotate charges in the Martini resolution part if there are any
- draw your string to check it's correctness

Consider the difference between the Martini 3 mappings for Toluene and
2,4-dichlorotoluene shown below.

.. list-table:: Martini Mappings of Toluene and 2,4-dichlorotoluene
   :widths: auto
   :header-rows: 0

   * - .. image:: /images/toluene.jpeg
         :width: 400px
     - .. image:: /images/dichlorotoluene.jpeg
         :width: 400px

In Toluene the two TC5 beads are equivalent and connect the same way to the SC4
bead. Therefore the CGsmiles string below is valid.

.. code::

   Toluene
   {[#SC4]1[#TC5][#TC5]1}.{#SC4=Cc(c[!])c[!],#TC5=[!]ccc[!]}

However, in 2,4-dichlorotoluene there are two SX3 beads which are equivalent
except for the fact that they connect to the SC4 beads at different carbons. Once
the carbon with the chlorine connects to the SC4 and once the carbon without the
chlorine connects. To represent this connectivity correctly you need to use two
different fragments and two differently labeled bonding operators in the CGsmiles 
string as shown below:

.. code::

   2,4-dichlorotoluene
   {[#SC4]1[#SX3][#SX3A]1}.{#SC4=Cc[$a]c[$],#SX3=Clc[$a]c[$b],#SX3A=Clc[$b]c[$]}


Examples
--------
An extensive list with examples can be found in the publication.
