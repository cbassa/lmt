#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <getopt.h>
#include "filterbank.h"

#define LIM 256

int main(int argc,char *argv[])
{
  int arg=0;
  int i=0,j,k=0,nsub=128;
  fftwf_complex *rp1,*rp2;
  FILE *infile,*outfile;
  char infname[LIM],outfname[LIM];
  int bytes_read;
  float *ap;
  struct filterbank fbin,fbout;

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
      nsub=(unsigned int) atoi(optarg);
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
  bytes_read=fread(&fbin,1,sizeof(struct filterbank),infile);

  // Copy filterbank struct
  fbout.mjd_start=fbin.mjd_start;
  fbout.intmjd=fbin.intmjd;
  fbout.intsec=fbin.intsec;
  strcpy(fbout.source,fbin.source);
  strcpy(fbout.telescope,fbin.telescope);
  strcpy(fbout.instrument,fbin.instrument);
  fbout.freq=fbin.freq;
  fbout.bw=fbin.bw;
  fbout.npol=1; // Single polarization
  fbout.nbit=fbin.nbit; // Floats
  fbout.ndim=1; // Real output
  fbout.nchan=fbin.nchan; // Number of channels
  fbout.fsamp=fbin.fsamp; // Channelsize
  fbout.tsamp=fbin.tsamp*nsub; // Updated sample size

  // Write header struct
  fwrite(&fbout,1,sizeof(struct filterbank),outfile);

  // Allocate buffers
  rp1=fftwf_malloc(sizeof(fftwf_complex)*fbin.nchan);
  rp2=fftwf_malloc(sizeof(fftwf_complex)*fbin.nchan);
  ap=(float *) malloc(sizeof(float)*fbin.nchan);

  // Read buffers
  for (;;) {
    // Reset buffers
    for (j=0;j<fbin.nchan;j++) 
      ap[j]=0.0;

    // Accumulate subints
    for (i=0;i<nsub;i++) {
      bytes_read=fread(rp1,sizeof(fftwf_complex),fbin.nchan,infile);
      bytes_read=fread(rp2,sizeof(fftwf_complex),fbin.nchan,infile);

      // Sum results
      for (j=0;j<fbin.nchan;j++) {
	ap[j]+=rp1[j][0]*rp1[j][0]+rp1[j][1]*rp1[j][1];
	ap[j]+=rp2[j][0]*rp2[j][0]+rp2[j][1]*rp2[j][1];
      }
    }

    // Exit when buffer is empty
    if (bytes_read==0)
      break;

    // Scale
    for (j=0;j<fbin.nchan;j++) 
      ap[j]/=(float) nsub;

    // Write
    fwrite(ap,sizeof(float),fbin.nchan,outfile);
  } 

  // Close
  fclose(infile);
  fclose(outfile);

  // Free
  fftwf_free(rp1);
  fftwf_free(rp2);
  free(ap);

  return 0;
}
