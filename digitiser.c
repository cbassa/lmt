// Repack the complex voltige after FFT into time-domain

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <stdint.h>
#include "timeseries.h"

#define LIM 256

int main(int argc,char *argv[])
{
  int arg=0;
  FILE *infile,*outfile,*dadahdr;
  char infname[LIM],outfname[LIM],hdrfname[LIM];
  char *buffer,ut[30];
  struct timeseries ts;
  unsigned int bytes_read,blocksize=64000;

  // Decode options
  while ((arg=getopt(argc,argv,"i:o:b:"))!=-1) {
    switch (arg) {

    case 'i':
      strcpy(infname,optarg);
      break;

    case 'o':
      strcpy(outfname,optarg);
      break;

    case 'b':
      blocksize=(unsigned int) atoi(optarg);
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

  // Read in block header
  fread(&ts,1,sizeof(struct timeseries),infile);

  // Open output file
  outfile=fopen(outfname,"w");

  // Check if output file exists
  if (outfile==NULL) {
    fprintf(stderr,"Error opening %s\n",outfname);
    exit;
  }

  // Write the header
  fwrite(&ts,1,sizeof(struct timeseries),outfile);

  /**** Loop over? ****/

  // Read the data
  buffer=(*char) malloc(sizeof(char)*blocksize);

  // Repack the data, bit from ts.nbit
  for (i=0; i<N; i++) {
    if (Pol1[i]*scale>127) Pol1[i] = 127./scale;
    if (Pol2[i]*scale>127) Pol2[i] = 127./scale;
    if (Pol1[i]*scale<-128) Pol1[i] = -128./scale;
    if (Pol2[i]*scale<-128) Pol2[i] = -128./scale;
    Data[2*i] = (char) floor(Pol1[i]*scale + 0.);
    Data[2*i+1] = (char) floor(Pol2[i]*scale + 0.);
  }

  // Write the data

  /*************/

  // Clean up
  fclose(infile);
  fclose(outfile);
  free(buffer);

  return 0;
}
