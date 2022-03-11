export MFEM_INSTALL_PATH=/lore/hakimm2/opt/mfem_no_pumi
export HYPRE_DIR=$HYPRE_ROOT
export METIS_DIR=$METIS_ROOT
flags="-g -O0"
cmake .. \
  -DCMAKE_BUILD_TYPE=Debug \
  -DMFEM_PREFIX=$MFEM_INSTALL_PATH \
  -DUSE_OMEGAH=NO \
  -DCMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
  -DCMAKE_INSTALL_PREFIX:PATH="$CMAKE_BINARY_DIR/../bin" \
  -DCMAKE_C_COMPILER="mpicc" \
  -DCMAKE_CXX_COMPILER="mpicxx" \
  -DCMAKE_C_FLAGS="${flags}" \
  -DCMAKE_CXX_FLAGS="${flags}" \
  -DCMAKE_EXE_LINKER_FLAGS="-lpthread ${flags}"
