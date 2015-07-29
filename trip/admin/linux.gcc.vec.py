CC = 'gcc', 
CCFLAGS = '-O3 -mfpmath=sse -ffast-math -falign-functions=16 -fstrict-aliasing -ftree-vectorize -ftree-vectorizer-verbose=5 -fprefetch-loop-arrays --param vect-max-version-for-alias-checks=200 -pedantic -Wunused -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE'
CFLAGS = '-std=c99 -Wimplicit'
CXX = 'g++'
