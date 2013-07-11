CC = 'mpicc' 
CCFLAGS = '-O3 -xSSE4.2 -falign-functions=16 -restrict -vec-report3 -opt-report 3 -opt-report-phase=hlo -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -DIWAVE_USE_MPI'
CFLAGS = '-std=c99'
CXX = 'mpicxx'
