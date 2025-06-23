ARG MINIFORGE_VERSION=24.11.3-2
FROM condaforge/miniforge3:${MINIFORGE_VERSION} AS build
WORKDIR /cadet

USER root

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get -y install build-essential cmake libhdf5-dev libsuperlu-dev intel-mkl git git-lfs libeigen3-dev && \
    apt-get clean

COPY . CADET-Core
WORKDIR CADET-Core

RUN mkdir -p build

WORKDIR build

SHELL ["/bin/bash", "-c"]

ENV MKLROOT=/opt/intel/mkl

RUN cmake -DCMAKE_INSTALL_PREFIX="../install" -DENABLE_STATIC_LINK_DEPS=ON -DENABLE_STATIC_LINK_LAPACK=ON -DBLA_VENDOR=Intel10_64lp_seq ../

RUN make -j $(lscpu | grep 'CPU(s)' | head -n 1 | cut -d ':' -f2 | tr -d ' ') install

RUN /cadet/CADET-Core/install/bin/createLWE -o /cadet/CADET-Core/install/bin/LWE.h5
RUN /cadet/CADET-Core/install/bin/cadet-cli /cadet/CADET-Core/install/bin/LWE.h5

FROM condaforge/miniforge3:${MINIFORGE_VERSION} AS deploy
COPY --from=build /cadet/CADET-Core/install /cadet/CADET-Core/install
COPY --from=build /usr/lib/x86_64-linux-gnu/libsz.so.2 /cadet/CADET-Core/install/lib
ENV PATH="$PATH:/cadet/CADET-Core/install/bin"

RUN apt-get update && \
    apt-get -y install libhdf5-dev libsuperlu-dev git git-lfs && \
    apt-get clean

# Ensure CADET-Core still works
RUN cadet-cli --version

# CADET_PYTHON_HASH can be a branch name or a commit hash or vX.X.X
ARG CADET_PYTHON_HASH=master
RUN git clone https://github.com/cadet/CADET-Python /cadet/CADET-Python
WORKDIR /cadet/CADET-Python
RUN git checkout ${CADET_PYTHON_HASH}
RUN pip install .

# CADET_PROCESS_HASH can be a branch name or a commit hash or vX.X.X
ARG CADET_PROCESS_HASH=master
RUN git clone https://github.com/fau-advanced-separations/CADET-Process /cadet/CADET-Process
WORKDIR /cadet/CADET-Process
RUN git checkout ${CADET_PROCESS_HASH}
RUN pip install .

RUN pip install cadet-rdm

WORKDIR /tmp
