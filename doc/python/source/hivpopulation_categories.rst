hivpopulation Overview
======================
.. currentmodule:: FFPopSim

``hivpopulation`` is a subclass of ``haploid_highd``, hence inherits its methods.
Moreover, ``hivpopulation`` has some HIV-specific methods.

``hivpopulation`` has two phenotypic traits, replication and resistance.
They contribute linearly to the viral fitness. The relative weight is set by
the attribute ``treatment``.

.. note:: Clicking on the name of the function will take you to a more detailed explanation listing all arguments.

Initialization
^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   hivpopulation.__init__

Drug Treatment
^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   hivpopulation.treatment

Replication and Resistance
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   hivpopulation.set_trait_landscape
   hivpopulation.set_replication_landscape
   hivpopulation.set_resistance_landscape

   hivpopulation.get_replication_additive
   hivpopulation.get_resistance_additive
   hivpopulation.set_replication_additive
   hivpopulation.set_resistance_additive

I/O Convenience Functions
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   hivpopulation.read_replication_coefficients
   hivpopulation.read_resistance_coefficients
   hivpopulation.write_genotypes
   hivpopulation.write_genotypes_compressed
