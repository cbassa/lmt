#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <stdint.h>
#include "timeseries.h"
#include "lib/delays.h"

#define LIM 256

void Usage()
{
  printf("Usage: reader -i <input file> -o <output file> [-b <blocksize> (64000)] -c <calcfile> -g <gpsfile> -r <refgpsfile> -s <samples to skip>\n");
}

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
  long binnumber, reducedbinnumber;
  struct timeseries hdr;
  unsigned int bytes_read,blocksize=64000;
  int iset=0, oset=0, cset=0, gset=0, rset=0, sset=0;

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
      // Get number of samples to shift
    case 's':
      delays.skipsamples = atoi(optarg);
      sset = 1;
      break;
      
    default:
      return 0;
    }
  }
  
  // Check that all necessary parameters are set.
  if (!iset){ fprintf(stderr, "Please provide input parameter with -i\n"); Usage(); exit(0);}
  if (!oset){ fprintf(stderr, "Please provide output parameter with -o\n"); Usage(); exit(0);}
  if (!cset){ fprintf(stderr, "Please provide calcfilename with -c\n"); Usage(); exit(0);}
  if (!gset){ fprintf(stderr, "Please provide gps filename with -g\n"); Usage(); exit(0);}
  if (!rset){ fprintf(stderr, "Please provide reference gps filename with -r\n"); Usage(); exit(0);}
  if (!sset){ fprintf(stderr, "Please provide number of bins to shift with -s\n"); Usage(); exit(0);}

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
  if (hdr.ndim == 1) hdr.tsamp *= 2; // Sampling time of real data will become twice the value once complex.

  // Write header struct
  fwrite(&hdr,1,sizeof(struct timeseries),outfile);

  // Read the GPS clockfile
  ReadGPSFiles(&gps, hdr.mjd_start); 
  printf("clockoffset: %lf\n", gps.ClockOffset);

  // Read Calc for geometric delay correction
  calc.Poly = (double*) calloc(NPOLYS, sizeof(double));
  ReadCalcfile(&calc, hdr.mjd_start, delays.skipsamples, hdr.tsamp, 1);
  binnumber = (long)(calc.Offset/hdr.tsamp + 0.5); // binnumber from start of calc in output bins
  reducedbinnumber = binnumber-(long)(calc.CountCalcSteps*CALCSTEPSIZE/hdr.tsamp); // binnumber from start of current Calc step
  printf("calc offset: %e\n", calc.Offset);
  printf("binnumber: %ld\n", binnumber);
  printf("reducedbinnumber: %ld\n", reducedbinnumber);
  printf("samptime: %e\n", hdr.tsamp);
  startgeooffset = GetStartGeoOffset(&calc, hdr.tsamp, reducedbinnumber);
  printf("startgeooffset: %e s\n", startgeooffset);
  // Intialize delays that keep track on how many bins have been shifted
  delays.samplesshifted = delays.geoshifted = 0;
  Add_Geo_Delays(&delays, startgeooffset, hdr.tsamp);
  Add_Clock_Delays(&delays, &gps, hdr.tsamp);
  printf("skipsamples: %ld\n", delays.skipsamples);
  printf("samplesshifted: %ld\n", delays.samplesshifted);
  printf("geoshifted: %ld\n", delays.geoshifted);

  // Skipsamples in infile
  fseek(infile, 4*delays.skipsamples, SEEK_SET);

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
