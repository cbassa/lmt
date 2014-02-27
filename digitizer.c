#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <stdint.h>
#include "timeseries.h"

#define LIM 256

char digitize(float f,float scale,float offset)
{
  char c;

  // Apply scaling and offset
  f=f*scale+offset;

  // Keep within bounds
  if (f>127.0)
    f=127.0;
  else if (f<-128.0)
    f=-128.0;

  // Convert to 8bit char, adding to ensure correct rounding
  c=(char) floor(f);

  return c;
}

int main(int argc,char *argv[])
{
  int i,arg=0;
  FILE *infile,*outfile,*dadahdr;
  char infname[LIM],outfname[LIM],hdrfname[LIM];
  float *inbuf;
  char *outbuf,ut[30];
  struct timeseries ts;
  unsigned int bytes_read;
  float scale=1.0,offset=0.0;

  // Decode options
  while ((arg=getopt(argc,argv,"i:o:S:O:"))!=-1) {
    switch (arg) {

    case 'i':
      strcpy(infname,optarg);
      break;

    case 'o':
      strcpy(outfname,optarg);
      break;

    case 'S':
      scale=atof(optarg);
      break;

    case 'O':
      offset=atof(optarg);
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

  // Check if output file exists
  if (outfile==NULL) {
    fprintf(stderr,"Error opening %s\n",outfname);
    exit;
  }

  // Read header
  bytes_read=fread(&ts,1,sizeof(struct timeseries),infile);

  // Update information
  ts.nbit=8; // Back to 8bit chars

  // Write the header
  fwrite(&ts,1,sizeof(struct timeseries),outfile);

  // Allocate
  inbuf=(float *) malloc(sizeof(float)*ts.nsamp);
  outbuf=(char *) malloc(sizeof(char)*ts.nsamp);

    // Loop over contents
  for (;;) {
    // Read buffers
    bytes_read=fread(inbuf,sizeof(float),ts.nsamp,infile);

    // Exit when buffer is empty
    if (bytes_read==0)
      break;

    // Apply scaling
    for (i=0;i<ts.nsamp;i++) 
      outbuf[i]=digitize(inbuf[i],scale,offset);

    // Write
    fwrite(outbuf,sizeof(char),ts.nsamp,outfile);
  } 
  
  // Clean up
  fclose(infile);
  fclose(outfile);
  free(inbuf);
  free(outbuf);

  return 0;
}
