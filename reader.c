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
  delays_type delays;
  double startMJD, samptime, startgeooffset;
  long skipbins, binnumber, reducedbinnumber;
  struct timeseries hdr;
  unsigned int bytes_read,blocksize=64000;
  int iset=0, oset=0, bset=0, cset=0, gset=0, rset=0, sset=0;

  // Decode options
  while ((arg=getopt(argc,argv,"i:o:b:c:g:r:s:"))!=-1) {
    switch (arg) {
      // Get input filename (dada file)
    case 'i':
      strcpy(infname,optarg);
      iset = 1;
      break;
      // Get output filename (fifo)
    case 'o':
      strcpy(outfname,optarg);
      oset = 1;
      break;
      // Get the blocksize
    case 'b':
      blocksize=(unsigned int) atoi(optarg);
      bset = 1;
      break;
      // Get calc filename      
    case 'c':
      calc.filename = optarg;
      cset = 1;
      break;
      // Get GPS filename of telescope
    case 'g':
      gps.filename = optarg;
      gset = 1;
      break;
      // Get GPS filename of reference telescope
    case 'r':
      gps.filenameref = optarg;
      rset = 1;
      break;
      // Get number of bins to shift
    case 's':
      skipbins = atoi(optarg);
      if (skipbins < 0) {fprintf(stderr, "Error, binshift cannot be negative\n");exit(0);}
      sset = 1;
      break;
      
    default:
      return 0;
    }
  }
  
  // Check that all necessary parameters are set.
  if (!iset){ fprintf(stderr, "Please provide input parameter with -i\n"); exit(0);}
  if (!oset){ fprintf(stderr, "Please provide output parameter with -o\n"); exit(0);}
  if (!bset){ fprintf(stderr, "Please provide blocksize with -b\n"); exit(0);}
  if (!cset){ fprintf(stderr, "Please provide calcfilename with -c\n"); exit(0);}
  if (!gset){ fprintf(stderr, "Please provide gps filename with -g\n"); exit(0);}
  if (!rset){ fprintf(stderr, "Please provide reference gps filename with -r\n"); exit(0);}
  if (!sset){ fprintf(stderr, "Please provide number of bins to shift with -s\n"); exit(0);}

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
  binnumber = (long)(calc.Offset/hdr.tsamp + 0.5); // binnumber from start of calc in output bins
  reducedbinnumber = binnumber-(long)(calc.CountCalcSteps*CALCSTEPSIZE/hdr.tsamp); // binnumber from start of current Calc step
  printf("binnumber: %e\n", binnumber);
  printf("reducedbinnumber: %e\n", reducedbinnumber);
  printf("samptime: %e\n", hdr.tsamp);
  startgeooffset = GetStartGeoOffset(&calc, hdr.tsamp, reducedbinnumber);
  printf("startgeooffset: %e\n");
  Add_Geo_Delays(&delays, startgeooffset, hdr.tsamp);
  Add_Clock_Delays(&delays, &gps, hdr.tsamp);

  exit(0);
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
