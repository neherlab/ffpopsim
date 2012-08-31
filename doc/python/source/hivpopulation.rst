.. _hivpopulation:

.. currentmodule:: FFPopSim

hivpopulation Class Reference
=============================
.. autoclass:: hivpopulation
   :members: treatment,
             write_genotypes_compressed,
             set_trait_landscape,
             set_replication_landscape,
             set_resistance_landscape,
             get_replication_additive,
             get_resistance_additive,
             set_replication_additive,
             set_resistance_additive

   .. automethod:: hivpopulation.__init__(N=0, rng_seed=0, mutation_rate=3e-5, coinfection_rate=1e-2, crossover_rate=1e-3)
   .. automethod:: hivpopulation.write_genotypes(filename, sample_size, gt_label='', start=0, length=0)
   .. automethod:: hivpopulation.read_replication_coefficients(filename)
   .. automethod:: hivpopulation.read_resistance_coefficients(filename)
