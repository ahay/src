#ifndef _sf_files_h
#define _sf_files_h

#include "c99.h"
#include "file.h"

int sf_filedims (sf_file file, /*@out@*/ int *n);
int sf_filesize (sf_file file);
int sf_leftsize (sf_file file, int dim);
void sf_cp(sf_file in, sf_file out);
void sf_rm(const char* filename, bool force, bool verb, bool inquire);

#endif

/* 	$Id: files.h,v 1.2 2003/09/29 14:34:55 fomels Exp $	 */
