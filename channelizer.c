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
  char *buffer;
  struct timeseries ts;
  struct filterbank fb;
  float *bp1,*bp2;
  fftwf_complex *rp1,*rp2,*cp1,*cp2;
  fftwf_plan ftp1,ftp2;
  unsigned int bytes_read,nchan=128;

  // Decode options
  while ((arg=getopt(argc,argv,"i:o:n:"))!=-1) {
    switch (arg) {

    case 'i':
      strcpy(infname,optarg);
      break;

    case 'o':
      strcpy(outfname,optarg);
      break;

    case 'n':
      nchan=(unsigned int) atoi(optarg);
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

  // Copy timeseries struct information to filterbank struct
  fb.mjd_start=ts.mjd_start;
  fb.intmjd=ts.intmjd;
  fb.intsec=ts.intsec;
  strcpy(fb.source,ts.source);
  strcpy(fb.telescope,ts.telescope);
  strcpy(fb.instrument,ts.instrument);
  fb.freq=ts.freq;
  fb.bw=ts.bw;
  fb.npol=ts.npol;
  fb.nbit=-32; // Floats
  fb.ndim=2; // Complex output
  fb.nchan=nchan; // Number of channels
  fb.fsamp=(double) fb.bw/(double) fb.nchan;

  // Allocate 
  if (ts.ndim==1) {
    // real to complex FFTs
    bp1=(float *) malloc(2*nchan*sizeof(float));
    bp2=(float *) malloc(2*nchan*sizeof(float));
    cp1=fftwf_malloc(sizeof(fftwf_complex)*(nchan+1));
    cp2=fftwf_malloc(sizeof(fftwf_complex)*(nchan+1));

    // FFTW Plans
    ftp1=fftwf_plan_dft_r2c_1d(2*nchan,bp1,cp1,FFTW_MEASURE);
    ftp2=fftwf_plan_dft_r2c_1d(2*nchan,bp2,cp2,FFTW_MEASURE);

    // Update sampling time
    fb.tsamp=2*ts.tsamp*nchan;
  } else if (ts.ndim==2) {
    // complex to complex FFTs
    rp1=fftwf_malloc(sizeof(fftwf_complex)*nchan);
    rp2=fftwf_malloc(sizeof(fftwf_complex)*nchan);
    cp1=fftwf_malloc(sizeof(fftwf_complex)*nchan);
    cp2=fftwf_malloc(sizeof(fftwf_complex)*nchan);

    // FFTW Plans
    ftp1=fftwf_plan_dft_1d(nchan,rp1,cp1,FFTW_FORWARD,FFTW_MEASURE);
    ftp2=fftwf_plan_dft_1d(nchan,rp2,cp2,FFTW_FORWARD,FFTW_MEASURE);

    // Update sampling time
    fb.tsamp=ts.tsamp*nchan;
  }

  // Allocate character buffer
  buffer=(char *) malloc(sizeof(char)*4*nchan);

  // Write header struct
  fwrite(&fb,1,sizeof(struct filterbank),outfile);

  // Loop over contents
  do {
    // Read buffer
    bytes_read=fread(buffer,1,4*nchan,infile);

    // Repack character buffer into FFTW buffers
    if (ts.ndim==1) {
      for(j=0,k=0;j<4*nchan;j+=2,k++) {
	bp1[k]=(float) buffer[j];
	bp2[k]=(float) buffer[j+1];
      }
    } else if (ts.ndim==2) {
      for (j=0,k=0;j<4*nchan;j+=4,k++) {
	rp1[k][0]=(float) buffer[j];
	rp1[k][1]=(float) buffer[j+1];
	rp2[k][0]=(float) buffer[j+2];
	rp2[k][1]=(float) buffer[j+3];
      }
    }

    // Perform Fast Fourier Transform
    fftwf_execute(ftp1);
    fftwf_execute(ftp2);
    
    // Unpack
    for (j=0;j<nchan;j++) {
      if (ts.ndim==1) {
	k=j;
      } else if (ts.ndim==2) { // Swap spectrum halves
	if (j<nchan/2)
	  k=j+nchan/2;
	else
	  k=j-nchan/2;
      }
      rp1[k][0]=cp1[j][0]/(float) nchan;
      rp1[k][1]=cp1[j][1]/(float) nchan;
      rp2[k][0]=cp2[j][0]/(float) nchan;
      rp2[k][1]=cp2[j][1]/(float) nchan;
    }

    // Write
    fwrite(rp1,sizeof(fftwf_complex),nchan,outfile);
    fwrite(rp2,sizeof(fftwf_complex),nchan,outfile);
  } while (bytes_read!=0);

  // Close
  fclose(infile);
  fclose(outfile);

  // Free
  free(buffer);
  fftwf_free(rp1);
  fftwf_free(rp2);
  fftwf_free(cp1);
  fftwf_free(cp2);
  fftwf_destroy_plan(ftp1);
  fftwf_destroy_plan(ftp2);

  return 0;
}
