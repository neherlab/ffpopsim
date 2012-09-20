.. _Install:

Install
=======
.. warning:: If you have not browsed the :doc:`Requirements` section already, please
             do so, to make sure you have installed all the requirements for
             FFPopSim.

The installation of the FFPopSim Python module is done via the programs Make and
distutils_. Please refer also to the `INSTALL` file if the instructions below do
not satisfy your needs (or generate errors!).

Using the binaries
^^^^^^^^^^^^^^^^^^
The simplest way to install FFPopSim is using the binaries provided for 32 and
64 bit Linux and Mac OSX (10.6+) systems. You can choose either way:

#. copy manually the files from ``build/<your arch>`` into a folder included in
   your ``PYTHONPATH``, where ``<your arch>`` is your architecture, i.e. Linux
   or Mac and 32 or 64 bit, or

#. if you have distutils_, install it system-wide, calling the following command
   *as a superuser*:

.. code-block:: bash

   make python-install

Neither of these strategies will involve any building, hence you do not need
GSL_ nor BOOST_ to install the binaries.


If you want to install FFPopSim using distutils_, but into a different location
than the standard third-party Python packages directory, you can call directly
(if needed, as superuser):

.. code-block:: bash

   python2.7 setup.py install --skip-build --install-lib=<target_dir>

where `<target_dir>` is the target installation directory.


Building FFPopSim locally
^^^^^^^^^^^^^^^^^^^^^^^^^
To build the module locally, call in your shell:

.. code-block:: bash

   make python

This call creates the binaries for your system. You can follow the instructions
above from this point on.

Testing FFPopSim
^^^^^^^^^^^^^^^^
To test whether FFPopSim is installed correctly (and inserted into your
``PYTHONPATH``), you can open Python2.7 in a new shell and call::

   import FFPopSim

If you do not get any errors, the installation is successful. You can proceed to
the :doc:`First steps with FFPopSim` section.

Troubleshooting
^^^^^^^^^^^^^^^
In case of problems with the installation, please check in the :doc:`Requirements`
section that you have all necessary run-time packages.

Please consult the file ``INSTALL`` in the main package folder for further help.

.. _GSL: http://www.gnu.org/software/gsl/
.. _BOOST: http://www.boost.org/
.. _distutils: http://docs.python.org/library/distutils.html
