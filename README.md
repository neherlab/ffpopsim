[![Build Status](https://travis-ci.org/iosonofabio/ffpopsim.svg?branch=master)](https://travis-ci.org/iosonofabio/ffpopsim)
[![Documentation Status](https://readthedocs.org/projects/ffpopsim/badge/?version=latest)](https://readthedocs.org/projects/ffpopsim/?badge=latest)
[![License](https://img.shields.io/badge/license-GPL3-blue.svg)](http://www.gnu.org/copyleft/gpl.html)

* License:	GPL3
* Author:	Richard Neher, Fabio Zanini
* Creation:	2012/09/12

Description
------------
![Genetic drift and draft](/data/drift_vs_draft.png)

**FFPopSim** is a **C++** and **Python** library to simulate large populations that are polymorphic at many loci. It allows for asexual and recombining populations and complex fitness functions, including pairwise and higher order epistasis. It is designed to study the effects of linked selection, the rare processes in large populations, and can be used to address a large variety of population genetics problems.

FFPopSim is easy to use via the Python interface and can be extended or modified both at the C++ or Python level.

Downloads
---------
The official sources and releases are online at:

https://github.com/neherlab/ffpopsim

https://github.com/iosonofabio/ffpopsim

Docs
----
Documentation and examples are online at:

http://webdav.tuebingen.mpg.de/ffpopsim/

http://ffpopsim.readthedocs.io/en/latest/

(The latter only for the Python interface.)

Install
-------

To install using pip:

```bash
pip install FFPopSim
```

To install the library, read the INSTALL file.
The two configuration files are:

- Makefile
- setup.py (for the Python extension)

Enjoy.
