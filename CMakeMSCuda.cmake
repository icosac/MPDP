SET(CUDA_ON "ON")
SET(CUDA_FLAGS "-O2 -std=c++11 -DCUDA_ON --default-stream per-thread --compiler-options -Wall -Wno-reorder")
SET(CUDA_ARCH "72;75")
SET(CUDA_SEP_COMP "ON")