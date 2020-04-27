#include <rsf.h>
#include <time.h>

#ifndef _BENCH_H

typedef struct bench_t bench_t;
/*^*/

struct bench_t{
  char const *fname;
  clock_t start;
  clock_t end;
  float total;
};
/*^*/

#endif

bench_t *timer[100];
int counter=0;

void tic(char const *name)
/*< start the timing >*/
{
  timer[counter]=malloc(sizeof(bench_t));
  timer[counter]->fname=name;
  timer[counter]->start=clock();
  counter++;
}

void toc(char const *name)
/*< end the timing >*/
{
  int ic=0;
  while (strcmp(timer[ic]->fname,name))
    ic++;
  timer[ic]->end=clock();
  timer[ic]->total=(float) (timer[ic]->end-timer[ic]->start) / CLOCKS_PER_SEC;
}

void printprof()
/*< print the info and clear the memory>*/
{
  char info[100];
  sf_warning("=========================================================== ");
  sf_warning("PROFILING: [CPU time] ");
  sf_warning("=========================================================== ");
  for (int ic=0; ic<counter; ic++){
    int len=strlen(timer[ic]->fname);
    int spacelen=50-len;
    sprintf(info,"%s",timer[ic]->fname);
    for (int is=0; is<spacelen; is++)
      sprintf(info,"%s ",info);
    sprintf(info,"%s: %7.3g [s]", info, timer[ic]->total);
    sf_warning(info);
    free(timer[ic]);
  }
  sf_warning("=========================================================== ");
  sf_warning("=========================================================== ");

}
