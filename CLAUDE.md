# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**FFPopSim** is a C++ library with Python bindings for simulating large populations that are polymorphic at many loci. It supports asexual and recombining populations with complex fitness functions, including pairwise and higher-order epistasis. The library is designed to study linked selection and rare processes in large populations.

## Build System

### C++ Library

Build the static library:
```bash
make src
```

Build and run C++ tests:
```bash
make tests      # Compile tests
make run-tests  # Execute all test binaries
```

The library depends on:
- GSL (GNU Scientific Library)
- Boost
- C++ compiler with O2 optimization by default

### Python Bindings

Build Python extension module:
```bash
make python
```

Install Python package system-wide (requires sudo):
```bash
make python-install
```

Run Python examples:
```bash
./run-examples                # Run all examples
./run-examples --skip-slow    # Skip slow examples
```

The Python bindings use SWIG to wrap the C++ library. Pre-generated SWIG wrappers are in `src/python/FFPopSim_wrap.cpp`. To regenerate SWIG files (rarely needed):
```bash
make swig
```

### Documentation

Build C++ documentation (requires Doxygen):
```bash
make doc
```

Build Python documentation (requires Sphinx):
```bash
make python-doc
# Or directly in doc/python:
cd doc/python && make html
```

Build documentation using Docker:
```bash
make docker-docs
```

### Full Build

Build everything (C++, Python, tests):
```bash
make
```

Build everything including documentation:
```bash
make all
```

### Clean Targets

```bash
make clean           # Clean built objects and binaries
make clean-all       # Clean everything including documentation
make clean-src       # Clean C++ library only
make clean-tests     # Clean test binaries only
make clean-python    # Clean Python builds
```

## Architecture

### Core Components

The library is organized into three main layers:

**1. Generic Layer** (`ffpopsim_generic.h`, `sample.cpp`)
- Common data structures and utilities
- `sample` class for statistical analysis
- Random number generation using GSL (gsl_rng_taus2)
- Recombination mode constants: `FREE_RECOMBINATION`, `CROSSOVERS`, `SINGLE_CROSSOVER`

**2. Low-Dimensional Module** (`ffpopsim_lowd.h`, `haploid_lowd.cpp`, `hypercube_lowd.cpp`)
- Simulations for small numbers of loci (L ≤ 20)
- `hypercube_lowd`: Represents fitness landscapes or genotype frequencies on binary hypercube
- `haploid_lowd`: Population evolution on the hypercube
- FFT-based transformations between function and coefficient representations
- Scales as O(2^L) for genotypes and O(3^L) for recombination intermediates

**3. High-Dimensional Module** (`ffpopsim_highd.h`, `haploid_highd.cpp`, `hypercube_highd.cpp`)
- Simulations for large numbers of loci (hundreds to thousands)
- `hypercube_highd`: Sparse representation using bitsets (boost::dynamic_bitset)
- `haploid_highd`: Population evolution with many loci
- `multiLocusGenealogy` and `rootedTree`: Genealogical tracking
- Uses word-based optimization with WORDLENGTH=28

**4. HIV-Specific Module** (`hivpopulation.h`, `hivpopulation.cpp`, `hivgene.cpp`)
- Specialized implementation for HIV evolution modeling
- HIV gene-specific structures and fitness models

### Key Design Patterns

- **Two Representations**: Low-dimensional hypercubes use both function values (genotype → fitness) and Fourier coefficients (epistatic effects). FFT transforms between them.
- **Sparse Storage**: High-dimensional module only stores observed genotypes and low-order epistatic coefficients.
- **State Management**: `hypercube_lowd` tracks whether `func` or `coeff` is up-to-date via `state` field (HC_FUNC, HC_COEFF, HC_FUNC_EQ_COEFF).

### Python Interface

Python bindings expose C++ classes directly:
- `FFPopSim.haploid_lowd(L)` - Low-dimensional population
- `FFPopSim.haploid_highd(L)` - High-dimensional population
- `FFPopSim.hivpopulation()` - HIV-specific population

Example workflow:
```python
import FFPopSim as h
pop = h.haploid_lowd(4)
pop.set_allele_frequencies([0, 0.3, 0.6, 0.9], N=1000)
pop.evolve(100)
pop.plot_diversity_histogram()
```

## Testing

### C++ Tests

Test binaries in `tests/`:
- `lowd` - Low-dimensional simulations
- `highd` - High-dimensional simulations
- `recombination_lowd` - Recombination in low-d
- `hivpopulation` - HIV-specific tests
- `test_genealogy` - Genealogical tracking

Run individual test:
```bash
./tests/lowd
```

### Python Tests

Python test examples are in `tests/python_*.py`. These are not automated unit tests but demonstration scripts.

## CI/CD

GitHub Actions workflow (`.github/workflows/build_and_publish.yml`):
1. Runs C++ tests on Ubuntu 22.04
2. Builds Python wheels for Linux/macOS and Python 3.6-3.11 using cibuildwheel
3. Creates source distribution (sdist)
4. Publishes to PyPI on tagged releases

## Development Environment

Docker development container:
```bash
./docker-dev  # Enter development environment with all dependencies
```

The Dockerfile installs:
- Debian 11.5 base
- Build tools (gcc, make)
- GSL and Boost libraries
- Python 3 and pip
- Development headers

## Configuration Files

- `Makefile` - Main build configuration (edit PYTHON, CXX, OPTIMIZATION_LEVEL at top)
- `setup.py` - Python extension build (edit include_dirs and library_dirs for GSL/Boost)
- `pyproject.toml` - Python package metadata
- `cibuildwheel.toml` - Wheel building configuration
- `.readthedocs.yml` - Read the Docs configuration

## Output Locations

All build artifacts go to `pkg/`:
- `pkg/lib/libFFPopSim.a` - Static C++ library
- `pkg/include/` - C++ headers
- `pkg/python/` - Python module (FFPopSim.py and _FFPopSim.so)
- `pkg/doc/cpp/` - C++ documentation
- `pkg/doc/python/` - Python documentation

## Common Issues

**GSL/Boost not found**: Edit `CXXFLAGS` in Makefile to add include paths (e.g., `-I/opt/local/include`) and edit `setup.py` to add library paths.

**Import errors**: Ensure `pkg/python` is in PYTHONPATH or install system-wide with `make python-install`.

**SWIG version**: Only SWIG 1.3 or 2.0 are officially supported (not 3.0+), but pre-generated wrappers are included.
