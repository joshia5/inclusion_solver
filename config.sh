export MFEM_INSTALL_PATH=path/to/mfem/install
export HYPRE_DIR=path/to/hypre/install
export METIS_DIR=path/to/metis/install
flags="-g -O0"
cmake .. \
  -DCMAKE_BUILD_TYPE=Debug \
  -DMFEM_PREFIX=$MFEM_INSTALL_PATH \
  -DUSE_OMEGAH=NO \
  -DCMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
  -DCMAKE_INSTALL_PREFIX:PATH="$CMAKE_BINARY_DIR/../bin" \
  -DCMAKE_C_COMPILER="/usr/bin/mpicc.mpich" \
  -DCMAKE_CXX_COMPILER="/usr/bin/mpicxx.mpich" \
  -DCMAKE_C_FLAGS="${flags}" \
  -DCMAKE_CXX_FLAGS="${flags}" \
  -DCMAKE_EXE_LINKER_FLAGS="-lpthread ${flags}"
