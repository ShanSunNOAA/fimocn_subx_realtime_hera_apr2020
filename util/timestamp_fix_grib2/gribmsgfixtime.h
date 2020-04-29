#ifndef TRAVERSEGRIB2_H
#define TRAVERSEGRIB2_H 1

#include <stdio.h>

long int gribmsgfixtime (FILE *fp);
long int checkgribmsgtime (FILE *fp);
void set_hrbytes (uint64_t fn_hr);

#endif	/* TRAVERSEGRIB2_H */
