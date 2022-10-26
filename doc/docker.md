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

- runs `docker build`, producing `neherlab/ffpopsim_builder` docker image, including all necessary dependencies and
  tools
- runs `docker run` command, starting a container based on `neherlab/ffpopsim_builder` image
- the arguments to the script are passed through to `docker run` (more specifically, the bash shell inside the
  container)
- the project directory is mounted into `/workdir` inside the container

This will run `make` inside the container, which builds both C++ and Python artifacts:

```bash
./docker-dev make
```

You can run arbitrary commands:

```bash
./docker-dev bash -c "PYTHONPATH='pkg/python' python -c 'import FFPopSim as h; pop = h.haploid_lowd(5); print pop'"
```

and arbitrary scripts:

```bash
./docker-dev python my_script.py
```

or start a bash shell into the container:

```bash
./docker-dev
```

For more commands, see the original dev documentation in the `INSTALL` file and in the docs website.

Note: the first build might take some time

Note: that these containers are throwaway and nothing persists inside them across runs. But you can persist files
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
  --user=1000:1000 \
  --volume=$(pwd):/workdir \
  --workdir=/workdir/ \
  --env=UID=1000 \
  --env=GID=1000 \
  --env=USER=user \
  --env=GROUP=user \
  neherlab/ffpopsim_builder \
  bash
```

Substitute your command instead of `bash` in the end.

