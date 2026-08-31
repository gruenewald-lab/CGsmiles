Installation
============

You can install **cgsmiles** from PyPI using pip:

.. code:: bash

   pip install cgsmiles

Alternatively, install **cgsmiles** from conda-forge using conda.
Make sure to check that the version is the most recent one:

.. code:: bash

   conda install -c conda-forge cgsmiles

To install the development version directly off GitHub use:

.. code:: bash

   pip install git+https://github.com/gruenewald-lab/CGsmiles.git

Note that some modules depend on optional dependencies. In particular,
the drawing module depends on the `scipy <https://scipy.org>`__
and `matplotlib <https://matplotlib.org>`__ packages. The rdkit
module depends on `rdkit <https://www.rdkit.org>`__.
These dependencies need to be installed before the module can be used.
