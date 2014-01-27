#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "filterbank.h"
#include <fftw3.h>
#include <getopt.h>

#define LIM 256

int main(int argc,char *argv[])
{
  int arg=0;
  int i=0,j,k=0,nchan,nsub=128;
  fftwf_complex *rp1,*rp2;
  FILE *infile,*outfile;
  char infname[LIM],outfname[LIM];
  int bytes_read;
  float *ap1,*ap2;
  struct filterbank fb;

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
  bytes_read=fread(&fb,1,sizeof(struct filterbank),infile);

  // Allocate buffers
  rp1=fftwf_malloc(sizeof(fftwf_complex)*fb.nchan);
  rp2=fftwf_malloc(sizeof(fftwf_complex)*fb.nchan);
  ap1=(float *) malloc(sizeof(float)*fb.nchan);
  ap2=(float *) malloc(sizeof(float)*fb.nchan);

  // Read buffers
  do {
    // Reset buffers
    for (j=0;j<fb.nchan;j++) {
      ap1[j]=0.0;
      ap2[j]=0.0;
    }

    // Accumulate subints
    for (i=0;i<nsub;i++) {
      bytes_read=fread(rp1,sizeof(fftwf_complex),fb.nchan,infile);
      bytes_read=fread(rp2,sizeof(fftwf_complex),fb.nchan,infile);

      // Sum results
      for (j=0;j<fb.nchan;j++) {
	ap1[j]+=rp1[j][0]*rp1[j][0]+rp1[j][1]*rp1[j][1];
	ap2[j]+=rp2[j][0]*rp2[j][0]+rp2[j][1]*rp2[j][1];
      }
    }

    // Scale
    for (j=0;j<fb.nchan;j++) {
      ap1[j]/=(float) nsub;
      ap2[j]/=(float) nsub;
      printf("%d %d %f %f\n",k,j,ap1[j],ap2[j]);
    }
    printf("\n");
    k++;
  } while (bytes_read!=0);

  // Close
  fclose(infile);
  fclose(outfile);

  // Free
  fftwf_free(rp1);
  fftwf_free(rp2);
  free(ap1);
  free(ap2);

  return 0;
}
