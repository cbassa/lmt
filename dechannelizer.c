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
  float *bp1,*bp2;
  unsigned int bytes_read;
  double offset;
  float *buffer;
  int ndim=2; // Default is complex sampling

  // Decode options
  while ((arg=getopt(argc,argv,"i:o:r"))!=-1) {
    switch (arg) {

    case 'i':
      strcpy(infname,optarg);
      break;

    case 'o':
      strcpy(outfname,optarg);
      break;

    case 'r':
      ndim=1;
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
  ts.ndim=ndim; 
  ts.nchan=1; 
  ts.nsamp=4*fb.nchan;

  // Allocate 
  cp1=fftwf_malloc(sizeof(fftwf_complex)*fb.nchan);
  cp2=fftwf_malloc(sizeof(fftwf_complex)*fb.nchan);
  if (ts.ndim==1) {
    // Allocate for complex to real FFTs
    rp1=fftwf_malloc(sizeof(fftwf_complex)*(fb.nchan+1));
    rp2=fftwf_malloc(sizeof(fftwf_complex)*(fb.nchan+1));
    bp1=(float *) malloc(2*fb.nchan*sizeof(float));
    bp2=(float *) malloc(2*fb.nchan*sizeof(float));

    // FFTW Plans
    ftp1=fftwf_plan_dft_c2r_1d(2*fb.nchan,rp1,bp1,FFTW_MEASURE);
    ftp2=fftwf_plan_dft_c2r_1d(2*fb.nchan,rp2,bp2,FFTW_MEASURE);

    ts.tsamp=0.5*fb.tsamp/(double) fb.nchan; // Updated sample size
  } else if (ts.ndim==2) {
    // Allocate for complex to complex FFTs
    rp1=fftwf_malloc(sizeof(fftwf_complex)*fb.nchan);
    rp2=fftwf_malloc(sizeof(fftwf_complex)*fb.nchan);

    // FFTW Plans
    ftp1=fftwf_plan_dft_1d(fb.nchan,rp1,cp1,FFTW_BACKWARD,FFTW_MEASURE);
    ftp2=fftwf_plan_dft_1d(fb.nchan,rp2,cp2,FFTW_BACKWARD,FFTW_MEASURE);

    ts.tsamp=fb.tsamp/(double) fb.nchan; // Updated sample size
  }

  // Allocate output buffer
  buffer=(float *) malloc(sizeof(float)*ts.nsamp);
  
  // Print information
  printf("Dechannelizer: %d channel(s), %g MHz per channel, %g us %s sampling\n",ts.nchan,ts.bw,ts.tsamp*1e6,(ts.ndim==1 ? "real" : "complex"));

  // Write header struct
  fwrite(&ts,1,sizeof(struct timeseries),outfile);

  // Loop over contents
  for (;;) {
    // Read buffers
    bytes_read=fread(cp1,sizeof(fftwf_complex),fb.nchan,infile);
    bytes_read=fread(cp2,sizeof(fftwf_complex),fb.nchan,infile);

    // Exit when buffer is empty
    if (bytes_read==0)
      break;

    // Swap spectrum halves
    for (i=0;i<fb.nchan;i++) {
      if (ts.ndim==1) {
	j=i;
      } else if (ts.ndim==2) {
	if (i<fb.nchan/2)
	  j=i+fb.nchan/2;
	else
	  j=i-fb.nchan/2;
      }
      rp1[j][0]=cp1[i][0];
      rp1[j][1]=cp1[i][1];
      rp2[j][0]=cp2[i][0];
      rp2[j][1]=cp2[i][1];
    }

    // UNTESTED, THIS NEEDS TO BE CHECKED
    // Set N+1 channels to zero for c2r transform
    if (ts.ndim==1) {
      rp1[fb.nchan][0]=0.0;
      rp1[fb.nchan][1]=0.0;
      rp2[fb.nchan][0]=0.0;
      rp2[fb.nchan][1]=0.0;
    }

    // Perform Fast Fourier Transform
    fftwf_execute(ftp1);
    fftwf_execute(ftp2);
    
    // Repack into interleaving values
    if (ts.ndim==1) {
      for (i=0,j=0;i<2*fb.nchan;i++,j+=2) {
	buffer[j]=bp1[i];
	buffer[j+1]=bp2[i];
      }
    } else if (ts.ndim==2) {
      for (i=0,j=0;i<fb.nchan;i++,j+=4) {
	buffer[j]=cp1[i][0];
	buffer[j+1]=cp1[i][1];
	buffer[j+2]=cp2[i][0];
	buffer[j+3]=cp2[i][1];
      }
    }

    // Write
    fwrite(buffer,sizeof(float),ts.nsamp,outfile);
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
