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
  struct timeseries ts;
  unsigned int bytes_read,blocksize=64000;
  double delay=0.0;
  uint64_t samples_to_skip,bytes_to_skip;

  // Decode options
  while ((arg=getopt(argc,argv,"i:o:b:t:"))!=-1) {
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
      
    case 't':
      delay=atof(optarg);
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
  ts=read_dada_header(infile);
  // Set only parameter in ts that doesn't come from the header: the number of samples to put in the buffer and write out at once (this is reset by the dechannelizer, but setting it here means you can write out without channelizing and dechannelizing if you want to)
  ts.nsamp=blocksize;

  // Compute sample delay
  samples_to_skip=(uint64_t) (delay*ts.samples_per_second);
  bytes_to_skip=samples_to_skip*ts.bytes_per_sample;

  // Increment obs_offset
  ts.obs_offset+=bytes_to_skip;

  // Skip read
  fseek(infile,bytes_to_skip,SEEK_CUR);

  // Write header struct
  fwrite(&ts,1,sizeof(struct timeseries),outfile);

  // Print information
  printf("Reader: %s with %s at %s\n",ts.source,ts.instrument,ts.telescope);
  printf("Reader: %s timeseries, %g us sampling, %d polarizations, %d bits\n",(ts.ndim==1 ? "real" : "complex"),ts.tsamp*1e6,ts.npol,ts.nbit);
  printf("Reader: %g MHz bandwidth at %g MHz center frequency\n",ts.bw,ts.freq);
  if (delay>0.0)
    printf("Reader: skipping %lf seconds (%Ld samples, %Ld bytes)\n",delay,samples_to_skip,bytes_to_skip);

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
