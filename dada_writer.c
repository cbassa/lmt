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
  char infname[LIM],outfname[LIM],hdrfname[LIM]="Dada_header.txt";
  char *buffer,*header,ut[30];
  struct timeseries ts;
  unsigned int bytes_read;

  // Decode options
  while ((arg=getopt(argc,argv,"i:o:s:"))!=-1) {
    switch (arg) {

    case 'i':
      strcpy(infname,optarg);
      break;

    case 'o':
      strcpy(outfname,optarg);
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
  bytes_read=fread(&ts,1,sizeof(struct timeseries),infile);

  //Open header template
  dadahdr=fopen(hdrfname,"r");

  //Check if template file exists
  if (dadahdr==NULL) {
    fprintf(stderr,"Error opening header file %s\n",hdrfname);
    exit;
  }

  //Open output file
  outfile=fopen(outfname,"w");

  // Check if output file exists
  if (outfile==NULL) {
    fprintf(stderr,"Error opening %s\n",outfname);
    exit;
  }

  //Allocate memory
  header=(char *) malloc(sizeof(char)*HEADERSIZE);
  buffer=(char *) malloc(sizeof(char)*ts.nsamp);
  
  //Read in template
  bytes_read=fread(header,sizeof(char),HEADERSIZE,dadahdr);
  fclose(dadahdr);

  //Modify content based on info in struct
  ascii_header_set(header,"SOURCE","%s",ts.source);
  ascii_header_set(header,"TELESCOPE","%s",ts.telescope);
  ascii_header_set(header,"INSTRUMENT","%s",ts.instrument);
  ascii_header_set(header,"FREQ","%lf",ts.freq);
  ascii_header_set(header,"BW","%lf",ts.bw);
  ascii_header_set(header,"TSAMP","%lf",ts.tsamp*1.0e6);
  ascii_header_set(header,"FILE_SIZE","%ld",ts.file_size);
  ascii_header_set(header,"OBS_OFFSET","%ld",ts.obs_offset);
  ascii_header_set(header,"MJD_START","%.11lf",ts.mjd_start);
  mjd2date(ts.mjd_start,ut);
  ascii_header_set(header,"UTC_START","%s",ut);
  //ascii_header_set(header,"RA","%s",ra);
  //ascii_header_set(header,"DEC","%s",dec);
  //ascii_header_set(header,"FILE_NAME","%s",filebasename);

  /****************************************************/

  // Write header
  fwrite(header,sizeof(char),HEADERSIZE,outfile);

  // Iterate over file to read and write data block(s)
  for (;;) {
    // Read buffer
    bytes_read=fread(buffer,sizeof(char),ts.nsamp,infile);

    // Exit when buffer is empty
    if (bytes_read==0)
      break;

    // Write buffer
    fwrite(buffer,sizeof(char),bytes_read,outfile);
  } 

  //Clean up
  free(buffer);
  free(header);

  // Close files
  fclose(infile);
  fclose(outfile);

  return 0;
}
