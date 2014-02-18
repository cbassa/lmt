#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <getopt.h>
#include "timeseries.h"
#include "filterbank.h"

#define LIM 256

int main(int argc,char *argv[])
{
  int arg=0;
  int i,j,k;
  FILE *infile,*outfile;
  char infname[LIM],outfname[LIM];
  struct timeseries ts;
  struct filterbank fb;
  fftwf_complex *rp1,*rp2,*cp1,*cp2;
  fftwf_plan ftp1,ftp2;
  unsigned int bytes_read;
  double offset;

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
  outfile=fopen(outfname,"wb");
  
  // Check if output file exists
  if (outfile==NULL) {
    fprintf(stderr,"Error opening %s\n",outfname);
    exit;
  }

  // Read header
  bytes_read=fread(&fb,1,sizeof(struct filterbank),infile);

  // Copy filterbank struct to timeseries struct
  ts.mjd_start=fb.mjd_start;
  ts.intmjd=fb.intmjd;
  ts.intsec=fb.intsec;
  strcpy(ts.source,fb.source);
  strcpy(ts.telescope,fb.telescope);
  strcpy(ts.instrument,fb.instrument);
  ts.freq=fb.freq;
  ts.bw=fb.bw;
  ts.npol=fb.npol; 
  ts.nbit=fb.nbit; 
  ts.ndim=fb.ndim; 
  ts.nchan=1; 
  ts.tsamp=fb.tsamp/(double) fb.nchan; // Updated sample size

  // Allocate for complex to complex FFTs
  rp1=fftwf_malloc(sizeof(fftwf_complex)*fb.nchan);
  rp2=fftwf_malloc(sizeof(fftwf_complex)*fb.nchan);
  cp1=fftwf_malloc(sizeof(fftwf_complex)*fb.nchan);
  cp2=fftwf_malloc(sizeof(fftwf_complex)*fb.nchan);

  // FFTW Plans
  ftp1=fftwf_plan_dft_1d(fb.nchan,rp1,cp1,FFTW_BACKWARD,FFTW_MEASURE);
  ftp2=fftwf_plan_dft_1d(fb.nchan,rp2,cp2,FFTW_BACKWARD,FFTW_MEASURE);
  
  // Print information
  printf("Dechannelizer: %d channel(s), %g MHz per channel, %g us sampling\n",ts.nchan,ts.bw,ts.tsamp*1e6);

  // Write header struct
  fwrite(&ts,1,sizeof(struct timeseries),outfile);

  // Loop over contents
  for (;;) {
    // Read buffers
    bytes_read=fread(rp1,sizeof(fftwf_complex),fb.nchan,infile);
    bytes_read=fread(rp2,sizeof(fftwf_complex),fb.nchan,infile);

    // Exit when buffer is empty
    if (bytes_read==0)
      break;

    // Perform Fast Fourier Transform
    fftwf_execute(ftp1);
    fftwf_execute(ftp2);

    // Write
    fwrite(cp1,sizeof(fftwf_complex),fb.nchan,outfile);
    fwrite(cp2,sizeof(fftwf_complex),fb.nchan,outfile);
  } 

  // Close
  fclose(infile);
  fclose(outfile);

  // Free
  fftwf_free(rp1);
  fftwf_free(rp2);
  fftwf_free(cp1);
  fftwf_free(cp2);
  fftwf_destroy_plan(ftp1);
  fftwf_destroy_plan(ftp2);

  return 0;
}
