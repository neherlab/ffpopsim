# Build using Docker

### Install Docker

#### Linux

If you are on Linux, install "Docker Engine": https://docs.docker.com/engine/install/ubuntu/

It will run natively, without a VM.

A couple of additional handy commands:

Ubuntu: allow current user to run docker on Linux without sudo (adding the user into the `docker` group):

```bash
sudo groupadd docker
sudo gpasswd -a ${USER} docker

```

(reopen the shell for it to have effect)

Ubuntu: restart docker services:

```bash
sudo systemctl restart docker.{service,socket} containerd.service
```

#### Other OS

If you are on macOS or Windows, install "Docker Desktop": https://docs.docker.com/get-docker/

These OS don't have a Linux kernel, so they can only run Docker in a VM. Make sure to configure sufficient amount of RAM
and CPUs for the VM, in Docker Desktop settings.

### Build and run ffpopsim in Docker using dev script (recommended)

There is a convenience helper dev script `./docker-dev`, which:

- runs `docker build`, taking definition in the `Dockerfile` and producing docker image tagged `neherlab/ffpopsim_builder`, including all necessary dependencies and tools (the image is stores somewhere in the system directories and you only interact with it using `docker` commands and its tag).
- runs `docker run` command, starting a container based on `neherlab/ffpopsim_builder` image.
- the arguments to the script are passed through to `docker run` (more specifically, the bash shell inside the
  container).
- the project directory is mounted into `/workdir` inside the container.

Note: to run the script, you might need to install `readlink` utility. On Ubuntu: `apt-get install coreutils`

Run `make` inside the container, which builds both C++ and Python artifacts:

```bash
./docker-dev make
```

or if you want to use all CPU cores:

```bash
./docker-dev make -j $(nproc)
```

or to silence some of the warnings:

```
make -j $(nproc) CFLAGS='-w' CXXFLAGS='-w' SWIGFLAGS='-w302,511'
```

(but it is better to fix them when we have soem more time).


The results will be in `pkg/` dir. Notably, the python package is in the `pkg/python`. You can import it into python code by adding it to the `PYTHONPATH`, either in shell:

```
PYTHONPATH='pkg/python' python3 my_script.py
```

or in the caller script itself

```python
import sys
sys.path.append('pkg/python')
```


You can also use dev script to run arbitrary commands inside the container, including Python and bash:

```bash
./docker-dev bash -c "PYTHONPATH='pkg/python' python3 -c 'import FFPopSim as h; pop = h.haploid_lowd(5); print(pop)'"
```

and arbitrary scripts:

```bash
./docker-dev python3 my_script.py
```

a good starting point is examples in the examples directory:


```bash
./docker-dev python3 'examples/example.py'
```

You can also run all of them using this convenience script:

```bash
./docker-dev ./run-examples
```

(which is useful for testing your code changes)


You can start a long-running bash shell into the container:

```bash
./docker-dev bash
```

or start a long-running python shell into the container:

```bash
./docker-dev ipython
```


For more available commands, see the original dev documentation in the `INSTALL` file and in the docs website.

Note: the first run might take some time, because that's when the docker image is built. Same as if you modify the `Dockerfile` - the image is then rebuilt.

Note: these containers are throw-away and nothing persists inside them across runs. But you can persist files
to `/workdir` and they will appear in the project directory.

Note: docker access to the host filesystem is limited to the project directory, so you cannot go anywhere above it

### Build ffpopsim in Docker manually (if you feel adventurous or don't want to run dev script)

This section does the same as the previous, but without the bash script, and `docker` commands are spelled explicitly.

Build the docker image (once now, and every time you change Dockerfile) with:

```
docker build \
  --target=base \
  --tag=neherlab/ffpopsim_builder \ 
  --build-arg=UID=$(id -u) \
  --build-arg=GID=$(id -g) \
  --build-arg=USER=user \
  --build-arg=GROUP=user \
  .
```

And run it (as many times as you want) with:

```
docker run -it --rm --init \
  --name=neherlab-ffpopsim_builder \
  --user=$(id -u):$(id -u) \
  --volume=$(pwd):/workdir \
  --workdir=/workdir/ \
  --env=UID=$(id -u) \
  --env=GID=$(id -g) \
  --env=USER=user \
  --env=GROUP=user \
  neherlab/ffpopsim_builder \
  bash
```

Substitute your command instead of `bash` in the end.


### Run Jupyter Lab in docker

The dev script has a shortcut command `jupyter` or `j`, which starts a [Jupyter Lab](https://jupyter.org/) server from inside container, with access to all dependencies and to project directory:

```bash
./docker-dev j
```

Keep this command running and then open a browser at

```plaintext
http://localhost:15000
```

Run code interactively and draw plots:

```python
# Add FFPopSim to PYTHONPATH
import sys
sys.path.append('pkg/python')

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import FFPopSim as h

%matplotlib inline
%reload_ext autoreload
%autoreload 2

sns.set()

c = h.haploid_lowd(4)
c.set_allele_frequencies([0,0.3,0.6,0.9], N=1000) 
c.evolve(100)
c.plot_diversity_histogram()

```

You can install dependencies temporarily from the Lab cells:

```
import sys
!{sys.executable} -m pip install --quiet --upgrade numpy
```

Note, all Python code and shell commands in these notebooks are running inside the container. The container will be discarded as soon as it is stopped. Only changes to project directory will persist across runs, e.g. notebooks, checkpoints and generated data saved into project dir will stay. The project dir is mounted into container to the path `/workdir` and this is the default working directory.
