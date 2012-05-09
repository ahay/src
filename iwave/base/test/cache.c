/* courtesy of Max Deschantsreiter, April 2012
   exposes cache properties
*/

#include <unistd.h>
#include <stdio.h>
#include <math.h>

int main() {
  long size = sysconf(_SC_LEVEL1_ICACHE_SIZE) / pow(2, 10);
  long assoc = sysconf(_SC_LEVEL1_ICACHE_ASSOC);
  long line = sysconf(_SC_LEVEL1_ICACHE_LINESIZE);
  printf("level 1 icache size = %ldK, assoc = %ld, line size = %ld\n",
	 size, assoc, line);

  size = sysconf(_SC_LEVEL1_DCACHE_SIZE) / pow(2, 10);
  assoc = sysconf(_SC_LEVEL1_DCACHE_ASSOC);
  line = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);
  printf("level 1 dcache size = %ldK, assoc = %ld, line size = %ld\n",
	 size, assoc, line);

  size = sysconf(_SC_LEVEL2_CACHE_SIZE) / pow(2, 10);
  assoc = sysconf(_SC_LEVEL2_CACHE_ASSOC);
  line = sysconf(_SC_LEVEL2_CACHE_LINESIZE);
  printf("level 2 cache size = %ldK, assoc = %ld, line size = %ld\n",
	 size, assoc, line);

  return 0;
}
