.. _hivpopulation:

.. currentmodule:: FFPopSim

hivpopulation Class Reference
=============================
.. autoclass:: hivpopulation
   :members:

   .. automethod:: hivpopulation.__init__(N=0, rng_seed=0, mutation_rate=3e-5, coinfection_rate=1e-2, crossover_rate=1e-3)
   .. automethod:: hivpopulation.__str__
   .. automethod:: hivpopulation.__repr__
   .. autoattribute:: hivpopulation.treatment
   .. automethod:: hivpopulation.write_genotypes(filename, sample_size, gt_label='', start=0, length=0)
   .. automethod:: hivpopulation.write_genotypes_compressed
   .. automethod:: hivpopulation.read_replication_coefficients(filename)
   .. automethod:: hivpopulation.read_resistance_coefficients(filename)
   .. automethod:: hivpopulation.set_trait_landscape
   .. automethod:: hivpopulation.set_replication_landscape
   .. automethod:: hivpopulation.set_resistance_landscape
   .. automethod:: hivpopulation.get_replication_additive
   .. automethod:: hivpopulation.get_resistance_additive
   .. automethod:: hivpopulation.set_replication_additive
   .. automethod:: hivpopulation.set_resistance_additive
