#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <getopt.h>
#include "filterbank.h"

#define LIM 256

void send_string(FILE *,char *);
void send_int(FILE *,char *,int);
void send_double(FILE *,char *,double);
void send_float(FILE *,char *,float);

int main(int argc,char *argv[])
{
  int arg=0;
  FILE *infile,*outfile;
  char infname[LIM],outfname[LIM];
  struct filterbank fb;
  int bytes_read;
  double fch1,fchan;
  float *ap;

  // Decode options
  while ((arg=getopt(argc,argv,"i:o:"))!=-1) {
    switch (arg) {

    case 'i':
      strcpy(infname,optarg);
      break;

    case 'o':
      strcpy(outfname,optarg);
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
  bytes_read=fread(&fb,1,sizeof(struct filterbank),infile);

  // Test input
  if (fb.nbit!=-32) {
    fprintf(stderr,"%d bit output not supported!\n",fb.nbit);
    return -1;
  }
  if (fb.npol!=1) {
    fprintf(stderr,"%d pol output not supported!\n",fb.npol);
    return -1;
  }

  // Compute channel parameters
  fch1=fb.freq+0.5*fabs(fb.bw);
  fchan=-fabs(fb.bw)/(double) fb.nchan;

  // Write sigproc filterbank header
  send_string(outfile,"HEADER_START");
  send_string(outfile,"source_name");
  send_string(outfile,fb.source);
  send_int(outfile,"telescope_id",1);
  send_int(outfile,"machine_id",1);
  send_int(outfile,"data_type",1);
  send_double(outfile,"src_raj",0.0);
  send_double(outfile,"src_dej",0.0);
  send_double(outfile,"az_start",0.0);
  send_double(outfile,"za_start",0.0);
  send_double(outfile,"fch1",fch1);
  send_double(outfile,"foff",fchan);
  send_int(outfile,"nchans",fb.nchan);
  send_int(outfile,"nbits",32);
  send_double(outfile,"tstart",fb.mjd_start);
  send_double(outfile,"tsamp",fb.tsamp);
  send_int(outfile,"nifs",fb.npol);
  send_string(outfile,"HEADER_END");

  // Allocate
  ap=(float *) malloc(sizeof(float)*fb.nchan);

  // Loop over file
  for (;;) {
    // Read buffer
    bytes_read=fread(ap,sizeof(float),fb.nchan,infile);

    // Exit when buffer is emtpy
    if (bytes_read==0)
      break;

    // Write
    fwrite(ap,sizeof(float),fb.nchan,outfile);
  }

  // Close output file
  fclose(outfile);

  // Close input file
  fclose(infile);

  // Free
  free(ap);

  return 0;
}

// Send string
void send_string(FILE *file,char *string)
{
  int len;
  len=strlen(string);
  fwrite(&len,sizeof(int),1,file);
  fwrite(string, sizeof(char),len,file);
  
  return;
}

// Send float
void send_float(FILE *file,char *name,float floating_point)
{
  send_string(file,name);
  fwrite(&floating_point,sizeof(float),1,file);

  return;
}

// Send double
void send_double(FILE *file,char *name, double double_precision)
{
  send_string(file,name);
  fwrite(&double_precision,sizeof(double),1,file);

  return;
}

// Send int
void send_int(FILE *file,char *name, int integer)
{
  send_string(file,name);
  fwrite(&integer,sizeof(int),1,file);

  return;
}
