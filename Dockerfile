FROM debian:10 as base

SHELL ["bash", "-euxo", "pipefail", "-c"]

RUN set -euxo pipefail >/dev/null \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq --yes \
&& apt-get install -qq --no-install-recommends --yes \
  apt-transport-https \
  bash \
  build-essential \
  ca-certificates \
  curl \
  git \
  gnupg \
  gsl-bin \
  libboost-all-dev \
  libgsl-dev \
  lsb-release \
  python \
  python-numpy \
  python-pip \
  sudo \
  swig \
  time \
>/dev/null \
&& apt-get clean autoclean >/dev/null \
&& apt-get autoremove --yes >/dev/null \
&& rm -rf /var/lib/apt/lists/*

ARG USER=user
ARG GROUP=user
ARG UID
ARG GID

ENV USER=$USER
ENV GROUP=$GROUP
ENV UID=$UID
ENV GID=$GID
ENV TERM="xterm-256color"
ENV HOME="/home/${USER}"
ENV PATH="/usr/lib/llvm-${CLANG_VERSION}/bin:${HOME}/.local/bin:${PATH}"

# Install Python dependencies
RUN set -euxo pipefail >/dev/null \
&& curl -fsSL https://bootstrap.pypa.io/pip/2.7/get-pip.py | python

RUN set -euxo pipefail >/dev/null \
&& pip install --index-url=https://pypi.python.org/simple/ --user \
  biopython==1.76 \
  matplotlib \
  numpy \
  pandas

RUN set -euxo pipefail >/dev/null \
&& \
  if [ -z "$(getent group ${GID})" ]; then \
    groupadd --system --gid ${GID} ${GROUP}; \
  else \
    groupmod -n ${GROUP} $(getent group ${GID} | cut -d: -f1); \
  fi \
&& \
  if [ -z "$(getent passwd ${UID})" ]; then \
    useradd \
      --system \
      --create-home --home-dir ${HOME} \
      --shell /bin/bash \
      --gid ${GROUP} \
      --groups sudo \
      --uid ${UID} \
      ${USER}; \
  fi \
&& sed -i /etc/sudoers -re 's/^%sudo.*/%sudo ALL=(ALL:ALL) NOPASSWD: ALL/g' \
&& sed -i /etc/sudoers -re 's/^root.*/root ALL=(ALL:ALL) NOPASSWD: ALL/g' \
&& sed -i /etc/sudoers -re 's/^#includedir.*/## **Removed the include directive** ##"/g' \
&& echo "foo ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers \
&& touch ${HOME}/.hushlogin \
&& chown -R ${UID}:${GID} "${HOME}"

USER ${USER}

WORKDIR "/workdir"
