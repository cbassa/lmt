#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <getopt.h>
#include "filterbank.h"
#include <cpgplot.h>

#define LIM 256

int main(int argc,char *argv[])
{
  int arg=0;
  int i,j;
  FILE *infile;
  char infname[LIM];
  struct filterbank fb;
  int bytes_read;
  float *ap,apmin,apmax,f;

  // Decode options
  while ((arg=getopt(argc,argv,"i:"))!=-1) {
    switch (arg) {

    case 'i':
      strcpy(infname,optarg);
      break;

    default:
      return 0;

    }
  }

  // Open input file
  infile=fopen(infname,"r");
  
  // Check if input file exists
  if (infile==NULL) {
    fprintf(stderr,"Error opening %s\n",infname);
    exit;
  }

  // Read header
  bytes_read=fread(&fb,1,sizeof(struct filterbank),infile);

  // Test input
  if (fb.nbit!=-32) {
    fprintf(stderr,"%d bit output not supported!\n",fb.nbit);
    return -1;
  }
  if (fb.npol!=1) {
    fprintf(stderr,"%d pol output not supported!\n",fb.npol);
    return -1;
  }

  // Allocate
  ap=(float *) malloc(sizeof(float)*fb.nchan);

  // Initialize pgplot

  cpgopen("?");

  // Loop over file
  for (;;) {
    // Read buffer
    bytes_read=fread(ap,sizeof(float),fb.nchan,infile);

    // Exit when buffer is emtpy
    if (bytes_read==0)
      break;

    // Find extrema
    for (i=0;i<fb.nchan;i++) {
      if (i==0) {
	apmin=ap[i];
	apmax=ap[i];
      } else {
	if (ap[i]<apmin) apmin=ap[i];
	if (ap[i]>apmax) apmax=ap[i];
      }
    }

    cpgenv(fb.freq-0.5*fb.bw,fb.freq+0.5*fb.bw,apmin,apmax,0,0);
    cpglab("Frequency (MHz)","Flux (arbitrary)"," ");

    for (i=0;i<fb.nchan;i++) {
      f=fb.freq+fb.bw*(float) i/(float) (fb.nchan-1)-0.5*fb.bw;
      if (i==0)
	cpgmove(f,ap[i]);
      else
	cpgdraw(f,ap[i]);
    }
  }

  // End pgplot
  cpgend();


  // Close input file
  fclose(infile);

  // Free
  free(ap);

  return 0;
}

