ARG MINIFORGE_VERSION=24.11.3-2
FROM condaforge/miniforge3:${MINIFORGE_VERSION} AS build
WORKDIR /cadet

USER root

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get -y install build-essential cmake libhdf5-dev libsuperlu-dev intel-mkl git git-lfs libeigen3-dev && \
    apt-get clean

COPY . CADET
WORKDIR CADET

RUN mkdir -p build

WORKDIR build

SHELL ["/bin/bash", "-c"]

ENV MKLROOT=/opt/intel/mkl

RUN cmake -DCMAKE_INSTALL_PREFIX="../install" -DENABLE_STATIC_LINK_DEPS=ON -DENABLE_STATIC_LINK_LAPACK=ON -DBLA_VENDOR=Intel10_64lp_seq ../

RUN make -j $(lscpu | grep 'CPU(s)' | head -n 1 | cut -d ':' -f2 | tr -d ' ') install

RUN /cadet/CADET/install/bin/createLWE -o /cadet/CADET/install/bin/LWE.h5
RUN /cadet/CADET/install/bin/cadet-cli /cadet/CADET/install/bin/LWE.h5

FROM condaforge/miniforge3:${MINIFORGE_VERSION} AS deploy
COPY --from=build /cadet/CADET/install /cadet/CADET/install
COPY --from=build /usr/lib/x86_64-linux-gnu/libsz.so.2 /cadet/CADET/install/lib

RUN apt-get update && \
    apt-get -y install libhdf5-dev libsuperlu-dev git git-lfs && \
    apt-get clean

WORKDIR /tmp

RUN /cadet/CADET/install/bin/cadet-cli --version
