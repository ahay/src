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
  long ncalls;
};
/*^*/

#endif

bench_t *timer[100];
int counter=0;

void tic(char const *name)
/*< start the timing >*/
{
  if (counter==0){
    timer[counter]=calloc(1,sizeof(bench_t));
    timer[counter]->fname=name;
    timer[counter]->start=clock();
    timer[counter]->ncalls=1;
    counter++;
  }
  else{
    int ic=0;
    while (ic<counter){
      if (strcmp(timer[ic]->fname,name)==0) break;
      ic++;
    }

    if (ic==counter){
      timer[ic]=calloc(1,sizeof(bench_t));
      timer[counter]->fname=name;
      timer[ic]->start=clock();
      timer[ic]->ncalls=1;
      counter++;
    }
    else{
      timer[ic]->start=clock();
      timer[ic]->ncalls++;
    }
  }

}

void toc(char const *name)
/*< end the timing >*/
{
  int ic=0;
  while (strcmp(timer[ic]->fname,name))
    ic++;

  timer[ic]->end=clock();
  timer[ic]->total+=(float) (timer[ic]->end-timer[ic]->start) / CLOCKS_PER_SEC;
}

void printprof()
/*< print the info and clear the memory>*/
{
  char info[100];
  sf_warning("====================================================================");
  sf_warning("PROFILING: [CPU time] ");
  sf_warning("====================================================================");
  sf_warning("#call-----Function -------------------------------------Runtime-----");
  for (int ic=0; ic<counter; ic++){
    sprintf(info,"%4ld    %s", timer[ic]->ncalls, timer[ic]->fname);
    int len=strlen(info);
    int spacelen=50-len;
    for (int is=0; is<spacelen; is++)
      sprintf(info,"%s ",info);
    sprintf(info,"%s: %7.2f [s]", info, timer[ic]->total);
    sf_warning(info);
    free(timer[ic]);
  }
  sf_warning("====================================================================");
  sf_warning("====================================================================");

}
