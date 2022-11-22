FROM debian:11.5 as base

SHELL ["bash", "-euxo", "pipefail", "-c"]

# Install system dependencies
RUN set -euxo pipefail >/dev/null \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq --yes \
&& apt-get install -qq --no-install-recommends --yes \
  apt-transport-https \
  bash \
  bison \
  build-essential \
  ca-certificates \
  ccache \
  curl \
  git \
  gsl-bin \
  libboost-all-dev \
  libgsl-dev \
  libpcre2-dev \
  libpcre3-dev  \
  lsb-release \
  parallel \
  python3 \
  python3-dev \
  python3-pip \
  sudo \
  time \
>/dev/null \
&& apt-get clean autoclean >/dev/null \
&& apt-get autoremove --yes >/dev/null \
&& rm -rf /var/lib/apt/lists/*

# Install swig v3
RUN set -euxo pipefail >/dev/null \
&& git clone --recursive --branch="v4.1.0" "https://github.com/swig/swig" "/tmp/build/swig" \
&& cd "/tmp/build/swig" \
&& ./autogen.sh \
&& ./configure --prefix="/usr/local" \
&& make -j $(nproc) \
&& make install \
&& rm -rf "/tmp/build/swig"

# Setup a non-root user, group and home directory.
# This user is also configured to run `sudo` without password (inside container).
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
ENV SHELL="/bin/bash"

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


# Switch to non-root user
USER ${USER}


# Install pip dependencies from `requirements.txt`.
# If you want to add a dependency, add it to `requirements.txt` in the project root.
COPY /requirements.txt /
RUN set -euxo pipefail >/dev/null \
&& pip3 install --user -r /requirements.txt


# Import matplotlib the first time to build the font cache.
RUN set -euxo pipefail >/dev/null \
&& python3 -c "import matplotlib.pyplot" \


USER ${USER}

WORKDIR "/workdir"
