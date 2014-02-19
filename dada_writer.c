//With input of header struct and data block, write out in dada format

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <stdint.h>
#include "timeseries.h"
#include "ascii_header.c"
#include "mjd2date.h"

#define LIM 256
#define HEADERSIZE 4096

int main(int argc,char *argv[])
{
  int arg=0;
  FILE *infile,*outfile,*dadahdr;
  char infname[LIM],outfname[LIM],hdrfname[LIM];
  char *buffer,ut[30];
  struct timeseries ts;
  unsigned int bytes_read,blocksize=64000;

  // Decode options
  while ((arg=getopt(argc,argv,"i:o:b:s:"))!=-1) {
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

    case 's':
      strcpy(hdrfname,optarg);
      break;

    default:
      return 0;

    }
  }

  //Open input file
  infile=fopen(infname,"r");

  //Check if input file exists
  if (infile==NULL) {
    fprintf(stderr,"Error opening %s\n",infname);
    exit;
  }

  //Read in block header
  fread(&ts,1,sizeof(struct timeseries),infile);

  //Open header template
  dadahdr=fopen(hdrfname,"r");

  //Check if template file exists
  if (dadahdr==NULL) {
    fprintf(stderr,"Error opening %s\n",hdrfname);
    exit;
  }

  //Open output file
  outfile=fopen(outfname,"w");

  // Check if output file exists
  if (outfile==NULL) {
    fprintf(stderr,"Error opening %s\n",outfname);
    exit;
  }

  /*********Construct header for data file***********/

  //Allocate memo
  buffer=(char *) malloc(sizeof(char)*HEADERSIZE);

  //Initialise
  memset(buffer,0,HEADERSIZE);

  //Read in template
  fread(buffer,sizeof(char),HEADERSIZE,dadahdr);
  fclose(dadahdr);

  //Modify content based on info in struct
  ascii_header_set(buffer,"SOURCE","%s",ts.source);
  ascii_header_set(buffer,"TELESCOPE","%s",ts.telescope);
  ascii_header_set(buffer,"INSTRUMENT","%s",ts.instrument);
  ascii_header_set(buffer,"FREQ","%lf",ts.freq);
  ascii_header_set(buffer,"BW","%lf",ts.bw);
  ascii_header_set(buffer,"TSAMP","%lf",ts.tsamp*1.0e6);
  ascii_header_set(buffer,"FILE_SIZE","%ld",ts.file_size);
  ascii_header_set(buffer,"OBS_OFFSET","%ld",ts.obs_offset);
  ascii_header_set(buffer,"MJD_START","%.11lf",ts.mjd_start);
  mjd2date(ts.mjd_start,ut);
  ascii_header_set(buffer,"UTC_START","%s",ut);
  //ascii_header_set(buffer,"RA","%s",ra);
  //ascii_header_set(buffer,"DEC","%s",dec);
  //ascii_header_set(buffer,"FILE_NAME","%s",filebasename);

  /****************************************************/

  // Write header
  fwrite(buffer,sizeof(char),HEADERSIZE,outfile);

  // Free up buffer
  free(buffer);

  // Allocate memo for data block
  buffer=(char *) malloc(sizeof(char)*blocksize);

  // Iterate over file to read and write data block(s)
  for (;;) {
    // Read buffer
    bytes_read=fread(buffer,sizeof(char),blocksize,infile);

    // Exit when buffer is empty
    if (bytes_read==0)
      break;

    // Write buffer
    fwrite(buffer,sizeof(char),bytes_read,outfile);
  } 

  //Clean up
  free(buffer);
  fclose(infile);
  fclose(outfile);

  return 0;
}
