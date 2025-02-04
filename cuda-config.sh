export MFEM_INSTALL_PATH=/lore/joshia5/develop/build-mfem-cuda-omegah
flags="-g -O0"
cmake .. \
 -DCMAKE_BUILD_TYPE=Debug \
 -DMFEM_PREFIX=$MFEM_INSTALL_PATH \
 -DCMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
 -DCMAKE_INSTALL_PREFIX:PATH=$PWD/install \
 -DCMAKE_CXX_COMPILER=mpicxx \
 -DCMAKE_CXX_FLAGS="${flags}" \
 -DCMAKE_EXE_LINKER_FLAGS="-lpthread ${flags}" \
 -DINCLUSION_SOLVER_CUDA_ARCH=75 \
 -DUSE_OMEGAH=YES 
