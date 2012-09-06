Install
=======
.. warning:: If you have not browsed the :doc:`Requirements` section already, please
             do so, to make sure you have installed all the requirements for
             FFPopSim.

The installation of the FFPopSim Python module is done via the programs Make and
distutils_. To build the module locally, call in your shell:

.. code-block:: bash

   make python

This command creates the files ``FFPopSim.py`` and ``_FFPopSim.so`` in the
folder ``pkg/python``. Those file contain the whole FFPopSim package and must
always be located in the same folder.

To install FFPopSim on the whole system and include it in your ``PYTHONPATH``
(which is needed to import the library from arbitrary locations), call the
following command **as a superuser**:

.. code-block:: bash

   make python-install

If FFPopSim has been installed correctly, you can proceed to the :doc:`First
steps with FFPopSim` section.

Troubleshooting
^^^^^^^^^^^^^^^
Please consult the file ``INSTALL`` in the main package folder for help in case
of problems.
