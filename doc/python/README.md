# FFPopSim Python documentation

Readthedocs project for FFPopSim

## Building the docs with Docker (recommended)

Once you have [Docker](https://docs.docker.com/get-docker/) installed, run from the root of the project:

    make docker-docs

The HTML files will appear in `docs/build/html/` (for manual inspection) and served on `http://localhost:8000`. The package `sphinx-autobuild` will watch the files, rebuild the HTML and reload the page in the browser on changes. 


## Building the docs locally

Build dependencies are managed with [Conda](https://conda.io).

Enter the docs directory:

    cd doc/python

Install them into an isolated environment named `org.neherlab.ffpopsim.docs` with:

    conda env create

Enter the environment with:

    conda activate org.neherlab.ffpopsim.docs

You can now build the documentation with:

    make html

which invokes Sphinx to build static HTML pages in `doc/python/build/html/`.

On some platforms you can view them in the default browser by running:

    open build/html/index.html

or

    xdg-open build/html/index.html


Alternatively, you can run

    make autobuild

The HTML files will also appear in `doc/python/build/html/` (for manual inspection) and served on `http://localhost:8000`. The `sphinx-autobuild` will watch the files, rebuild the HTML and reload the page in the browser on changes.

Leave the environment with:

    conda deactivate
