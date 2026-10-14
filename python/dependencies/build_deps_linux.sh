#!/usr/bin/env bash
set -euxo pipefail

PREFIX=/opt/cadet_deps
NPROC=${NPROC:-$(getconf _NPROCESSORS_ONLN || echo 2)}

# Cache lives in repo mount inside manylinux container
CACHE_ROOT=/project/.ci_cache
CACHE_DIR="$CACHE_ROOT/cadet_deps_linux"
STAMP="$CACHE_DIR/.complete"

mkdir -p "$CACHE_ROOT"

# If cache exists, restore and skip building
if [ -f "$STAMP" ]; then
  echo "Restoring cached deps from $CACHE_DIR to $PREFIX"
  rm -rf "$PREFIX"
  mkdir -p "$(dirname "$PREFIX")"
  cp -a "$CACHE_DIR" "$PREFIX"
  echo "Cached deps restored."
  exit 0
fi

mkdir -p "$PREFIX"
mkdir -p /tmp/cadet-deps-src
cd /tmp/cadet-deps-src

yum -y install git make perl-core pkgconfig libatomic \
  openssl-devel zlib-devel bzip2-devel xz-devel \
  gcc-toolset-12-gcc gcc-toolset-12-gcc-c++

python -m pip install -U pip
python -m pip install -U cmake ninja patchelf

hash -r
cmake --version

export CMAKE_GENERATOR=Ninja
export PATH="$PREFIX/bin:$PATH"
export CMAKE_PREFIX_PATH="$PREFIX"
export PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig:$PREFIX/lib64/pkgconfig"
export LD_LIBRARY_PATH="$PREFIX/lib:$PREFIX/lib64:${LD_LIBRARY_PATH:-}"
export CFLAGS="-fPIC"
export CXXFLAGS="-fPIC"

# 1) OpenBLAS
if [ ! -d openblas ]; then
  git clone --depth 1 https://github.com/xianyi/OpenBLAS.git openblas
fi
cd openblas
make -j"$NPROC" TARGET=NEHALEM USE_OPENMP=0 NO_SHARED=0 > /tmp/openblas_build.log 2>&1
make PREFIX="$PREFIX" install > /tmp/openblas_install.log 2>&1
cd ..

# 2) Eigen
if [ ! -d eigen ]; then
  git clone --depth 1 --branch 3.4.0 https://gitlab.com/libeigen/eigen.git eigen
fi
cmake -S eigen -B eigen-build -DCMAKE_INSTALL_PREFIX="$PREFIX"
cmake --build eigen-build -j"$NPROC"
cmake --install eigen-build

# 3) oneTBB
if [ ! -d oneTBB ]; then
  git clone --depth 1 https://github.com/oneapi-src/oneTBB.git
fi
cmake -S oneTBB -B tbb-build \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DTBB_TEST=OFF \
  -DTBB_STRICT=OFF
cmake --build tbb-build -j"$NPROC"
cmake --install tbb-build

# 4) HDF5
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

# 5) SuiteSparse
if [ ! -d SuiteSparse ]; then
  git clone --depth 1 https://github.com/DrTimothyAldenDavis/SuiteSparse.git
fi
cmake -S SuiteSparse -B suitesparse-build \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DBUILD_SHARED_LIBS=ON \
  -DSUITESPARSE_ENABLE_PROJECTS="suitesparse_config;amd;btf;camd;ccolamd;colamd;cholmod;ldl;klu;umfpack;spqr" \
  -DBLA_VENDOR=OpenBLAS \
  -DCMAKE_PREFIX_PATH="$PREFIX"
cmake --build suitesparse-build -j"$NPROC"
cmake --install suitesparse-build

# 6) SuperLU
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

# Save into cache for next run
echo "Saving deps into cache $CACHE_DIR"
rm -rf "$CACHE_DIR"
mkdir -p "$CACHE_ROOT"
cp -a "$PREFIX" "$CACHE_DIR"
touch "$STAMP"
