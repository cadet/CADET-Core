#!/usr/bin/env bash
# Builds native dependencies for CADET-Core wheel builds (cibuildwheel / manylinux)
# Installs everything into /opt/cadet_deps

set -euo pipefail

PREFIX=/opt/cadet_deps
NPROC=${NPROC:-$(getconf _NPROCESSORS_ONLN || echo 2)}
mkdir -p "$PREFIX"
mkdir -p /tmp/cadet-deps-src
cd /tmp/cadet-deps-src

# Build tools
yum -y install git make cmake3 ninja-build perl-core pkgconfig patchelf \
  openssl-devel zlib-devel bzip2-devel xz-devel

# Use cmake3 as cmake on manylinux
if command -v cmake3 >/dev/null 2>&1; then
  ln -sf "$(command -v cmake3)" /usr/local/bin/cmake
fi

export PATH="$PREFIX/bin:$PATH"
export CMAKE_PREFIX_PATH="$PREFIX"
export PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig:$PREFIX/lib64/pkgconfig"
export LD_LIBRARY_PATH="$PREFIX/lib:$PREFIX/lib64:${LD_LIBRARY_PATH:-}"
export CFLAGS="-fPIC"
export CXXFLAGS="-fPIC"

# 1) OpenBLAS (BLAS/LAPACK provider)
if [ ! -d openblas ]; then
  git clone --depth 1 https://github.com/xianyi/OpenBLAS.git openblas
fi
cd openblas
make -j"$NPROC" DYNAMIC_ARCH=1 USE_OPENMP=0 NO_SHARED=0
make PREFIX="$PREFIX" install
cd ..

# 2) Eigen (header-only, install config files)
if [ ! -d eigen ]; then
  git clone --depth 1 --branch 3.4.0 https://gitlab.com/libeigen/eigen.git eigen
fi
cmake -S eigen -B eigen-build -DCMAKE_INSTALL_PREFIX="$PREFIX"
cmake --build eigen-build -j"$NPROC"
cmake --install eigen-build

# 3) oneTBB (for TBB::TBB)
if [ ! -d oneTBB ]; then
  git clone --depth 1 https://github.com/oneapi-src/oneTBB.git
fi
cmake -S oneTBB -B tbb-build \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DTBB_TEST=OFF \
  -DTBB_STRICT=OFF
cmake --build tbb-build -j"$NPROC"
cmake --install tbb-build

# 4) HDF5 (C library)
if [ ! -d hdf5 ]; then
  git clone --depth 1 --branch hdf5_1_14_4 https://github.com/HDFGroup/hdf5.git
fi
cmake -S hdf5 -B hdf5-build \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DBUILD_SHARED_LIBS=ON \
  -DHDF5_BUILD_TOOLS=OFF \
  -DHDF5_BUILD_EXAMPLES=OFF \
  -DHDF5_BUILD_HL_LIB=OFF \
  -DHDF5_BUILD_CPP_LIB=OFF \
  -DHDF5_ENABLE_Z_LIB_SUPPORT=ON \
  -DHDF5_ENABLE_SZIP_SUPPORT=OFF
cmake --build hdf5-build -j"$NPROC"
cmake --install hdf5-build

# 5) SuiteSparse (provides UMFPACK)
if [ ! -d SuiteSparse ]; then
  git clone --depth 1 https://github.com/DrTimothyAldenDavis/SuiteSparse.git
fi
cmake -S SuiteSparse -B suitesparse-build \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DBUILD_SHARED_LIBS=ON \
  -DSUITESPARSE_ENABLE_PROJECTS="amd;camd;colamd;ccolamd;cholmod;umfpack;spqr;metis;btf;klu;ldl;suitesparseconfig" \
  -DBLA_VENDOR=OpenBLAS \
  -DCMAKE_PREFIX_PATH="$PREFIX"
cmake --build suitesparse-build -j"$NPROC"
cmake --install suitesparse-build

# 6) SuperLU (optional, but you need it for 2D models)
if [ ! -d superlu ]; then
  git clone --depth 1 https://github.com/xiaoyeli/superlu.git
fi
cmake -S superlu -B superlu-build \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DBUILD_SHARED_LIBS=ON \
  -DTPL_ENABLE_LAPACKLIB=ON \
  -DBLA_VENDOR=OpenBLAS \
  -DCMAKE_PREFIX_PATH="$PREFIX"
cmake --build superlu-build -j"$NPROC"
cmake --install superlu-build

echo "Deps installed into $PREFIX"
