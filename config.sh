export MFEM_INSTALL_PATH=/home/morteza/Sandbox/install/mfem
export HYPRE_DIR=/home/morteza/Sandbox/install/hypre
export METIS_DIR=/home/morteza/Sandbox/install/metis
flags="-g -O0"
cmake /home/morteza/Sandbox/inclusion_solver \
  -DCMAKE_BUILD_TYPE=Debug \
  -DMFEM_PREFIX=$MFEM_INSTALL_PATH \
  -DCMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
  -DCMAKE_INSTALL_PREFIX:PATH="$CMAKE_BINARY_DIR/../bin" \
  -DCMAKE_C_COMPILER="/usr/bin/mpicc.mpich" \
  -DCMAKE_CXX_COMPILER="/usr/bin/mpicxx.mpich" \
  -DCMAKE_C_FLAGS="${flags}" \
  -DCMAKE_CXX_FLAGS="${flags}" \
  -DCMAKE_EXE_LINKER_FLAGS="-lpthread ${flags}"
