#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

struct filterbank {
  double mjd_start;
  unsigned int intmjd,intsec;
  char source[20],telescope[20],instrument[20];
  double freq,bw,tsamp,fsamp;
  unsigned int nbit,ndim,npol,nchan;
};

