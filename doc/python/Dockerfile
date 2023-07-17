FROM continuumio/miniconda3:23.3.1-0

SHELL ["bash", "-euxo", "pipefail", "-c"]

ARG DEBIAN_FRONTEND=noninteractive
ARG USER=user
ARG GROUP=user
ARG UID
ARG GID

ENV TERM="xterm-256color"
ENV HOME="/home/user"

# Install required packages if running Debian or Ubuntu
RUN set -euxo pipefail >/dev/null \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq --yes \
&& apt-get install -qq --no-install-recommends --yes \
  apt-transport-https \
  bash \
  bash-completion \
  ca-certificates \
  curl \
  git \
  gnupg \
  make \
  sudo \
>/dev/null \
&& apt-get clean autoclean >/dev/null \
&& apt-get autoremove --yes >/dev/null \
&& rm -rf /var/lib/apt/lists/*

# Make a user and group
RUN set -euxo pipefail >/dev/null \
&& \
  if [ -z "$(getent group ${GID})" ]; then \
    addgroup --system --gid ${GID} ${GROUP}; \
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
&& sed -i /etc/sudoers -re 's/^%sudo.*/%sudo ALL=(ALL:ALL) NOPASSWD:ALL/g' \
&& sed -i /etc/sudoers -re 's/^root.*/root ALL=(ALL:ALL) NOPASSWD:ALL/g' \
&& sed -i /etc/sudoers -re 's/^#includedir.*/## **Removed the include directive** ##"/g' \
&& echo "%sudo ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers \
&& echo "${USER} ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers \
&& touch ${HOME}/.hushlogin \
&& chown -R ${UID}:${GID} "${HOME}"

RUN set -euxo pipefail >/dev/null \
  && chown -R ${USER}:${GROUP} ${HOME}

COPY environment.yml ${HOME}/src/

WORKDIR ${HOME}/src

RUN set -euxo pipefail >/dev/null \
  && cd ${HOME}/src \
  && conda install anaconda-client -n base \
  && conda env create -f environment.yml --name org.neherlab.ffpopsim.docs

USER ${USER}

RUN set -euxo pipefail >/dev/null \
  && conda init bash \
  && echo "conda activate org.neherlab.ffpopsim.docs" >> ${HOME}/.bashrc

ENV PATH="/opt/conda/envs/org.neherlab.ffpopsim.docs/bin:${PATH}"

CMD bash -c "set -euxo pipefail \
  && source ${HOME}/.bashrc \
  && cd doc/python \
  && rm -rf build \
  && make autobuild \
  "
