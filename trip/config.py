CC = 'mpicc' 
CCFLAGS = '-O3 -pedantic -Wunused -Wno-long-long -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -DIWAVE_USE_MPI'
CFLAGS = '-std=c99 -Wimplicit'
CXX = 'mpicxx'
