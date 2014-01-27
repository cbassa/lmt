#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

struct timeseries {
  double mjd_start;
  unsigned int intmjd,intsec;
  char source[20],telescope[20],instrument[20];
  double freq,bw,tsamp;
  unsigned int nbit,ndim,npol,nchan;
  uint64_t obs_offset,samples_per_second,bytes_per_second,file_size,bytes_per_sample;
};
struct timeseries read_dada_header(FILE *file);
