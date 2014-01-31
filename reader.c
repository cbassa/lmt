#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <stdint.h>
#include "timeseries.h"
#include "lib/calc.h"

#define LIM 256

int main(int argc,char *argv[])
{
  int arg=0;
  FILE *infile,*outfile;
  char infname[LIM],outfname[LIM];
  char *buffer;
  calc_type calc;
  double startMJD, samptime;
  long skipbins;
  struct timeseries hdr;
  unsigned int bytes_read,blocksize=64000;

  // Decode options
  while ((arg=getopt(argc,argv,"i:o:b:c:"))!=-1) {
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
      
    case 'c':
      calc.filename = optarg;

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

  // Read Calc for geometric delay correction
  ReadCalcfile(calc, startMJD, skipbins, samptime, 1);

  // Allocate buffer
  buffer=(char *) malloc(sizeof(char)*blocksize);

  // Iterate over file
  do {
    // Read buffer
    bytes_read=fread(buffer,sizeof(char),blocksize,infile);

    // Write buffer
    fwrite(buffer,sizeof(char),bytes_read,outfile);
  } while (bytes_read!=0);

  // Close output file
  fclose(outfile);

  // Close input file
  fclose(infile);

  // Free buffer
  free(buffer);

  return 0;
}
