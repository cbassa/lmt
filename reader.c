#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <stdint.h>
#include "timeseries.h"
#include "lib/delays.h"

#define LIM 256

int main(int argc,char *argv[])
{
  int arg=0;
  FILE *infile,*outfile;
  char infname[LIM],outfname[LIM];
  char *buffer;
  calc_type calc;
  gps_type gps;
  double startMJD, samptime;
  long skipbins;
  struct timeseries hdr;
  unsigned int bytes_read,blocksize=64000;

  // Decode options
  while ((arg=getopt(argc,argv,"i:o:b:c:g:r:s:"))!=-1) {
    switch (arg) {
      // Get input filename (dada file)
    case 'i':
      strcpy(infname,optarg);
      break;
      // Get output filename (fifo)
    case 'o':
      strcpy(outfname,optarg);
      break;
      // Get the blocksize
    case 'b':
      blocksize=(unsigned int) atoi(optarg);
      break;
      // Get calc filename      
    case 'c':
      calc.filename = optarg;
      break;
      // Get GPS filename of telescope
    case 'g':
      gps.filename = optarg;
      break;
      // Get GPS filename of reference telescope
    case 'r':
      gps.filenameref = optarg;
      break;
      // Get number of bins to shift
    case 's':
      skipbins = atoi(optarg);
      if (skipbins < 0) {fprintf(stderr, "Error, binshift cannot be negative\n");exit(0);}
      break;
      
    default:
      return 0;

    }
  }

  // Open input file
  infile=fopen(infname,"r");
  
  // Check if input file exists
  if (infile==NULL) {
    fprintf(stderr,"Error opening inputfile %s\n",infname);
    exit;
  }

  // Open output file
  outfile=fopen(outfname,"w");
  
  // Check if output file can be opened
  if (outfile==NULL) {
    fprintf(stderr,"Error opening outputfile %s\n",outfname);
    exit;
  }

  // Read header
  hdr=read_dada_header(infile);

  // Write header struct
  fwrite(&hdr,1,sizeof(struct timeseries),outfile);

  // Read the GPS clockfile
  ReadGPSFiles(&gps, hdr.mjd_start); 

  // Read Calc for geometric delay correction
  calc.Poly = (double*) calloc(NPOLYS, sizeof(double));
  ReadCalcfile(&calc, hdr.mjd_start, skipbins, hdr.tsamp, 1);

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