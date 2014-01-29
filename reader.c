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
  FILE *infile,*outfile;
  char infname[LIM],outfname[LIM];
  char *buffer;
  struct timeseries hdr;
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

  // Open output file
  outfile=fopen(outfname,"w");
  
  // Check if input file exists
  if (outfile==NULL) {
    fprintf(stderr,"Error opening %s\n",outfname);
    exit;
  }

  // Read header
  hdr=read_dada_header(infile);

  // Write header struct
  fwrite(&hdr,1,sizeof(struct timeseries),outfile);

  // Print information
  printf("Reader: %s with %s at %s\n",hdr.source,hdr.instrument,hdr.telescope);
  printf("Reader: %s timeseries, %g us sampling, %d polarizations, %d bits\n",(hdr.ndim==1 ? "real" : "complex"),hdr.tsamp*1e6,hdr.npol,hdr.nbit);
  printf("Reader: %g MHz bandwidth at %g MHz center frequency\n",hdr.bw,hdr.freq);
  // Allocate buffer
  buffer=(char *) malloc(sizeof(char)*blocksize);

  // Iterate over file
  for (;;) {
    // Read buffer
    bytes_read=fread(buffer,sizeof(char),blocksize,infile);

    // Exit when buffer is empty
    if (bytes_read==0)
      break;

    // Write buffer
    fwrite(buffer,sizeof(char),bytes_read,outfile);
  } 

  // Close output file
  fclose(outfile);

  // Close input file
  fclose(infile);

  // Free buffer
  free(buffer);

  return 0;
}
