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
  int i=0,j,k=0;
  fftwf_complex *rp1,*rp2,*cp1,*cp2;
  FILE *infile1,*infile2,*outfile;
  char infname1[LIM],infname2[LIM],outfname[LIM];
  int bytes_read;
  struct filterbank fbin1,fbin2,fbout;

  // Decode options
  while ((arg=getopt(argc,argv,"i:I:o:"))!=-1) {
    switch (arg) {

    case 'i':
      strcpy(infname1,optarg);
      break;

    case 'I':
      strcpy(infname2,optarg);
      break;

    case 'o':
      strcpy(outfname,optarg);
      break;

    default:
      return 0;

    }
  }

  // Open input file
  infile1=fopen(infname1,"r");
  
  // Check if input file exists
  if (infile1==NULL) {
    fprintf(stderr,"Error opening %s\n",infname1);
    exit;
  }

  // Open input file
  infile2=fopen(infname2,"r");
  
  // Check if input file exists
  if (infile2==NULL) {
    fprintf(stderr,"Error opening %s\n",infname2);
    exit;
  }

  // Open output file
  outfile=fopen(outfname,"w");
  
  // Check if output file exists
  if (outfile==NULL) {
    fprintf(stderr,"Error opening %s\n",outfname);
    exit;
  }

  // Read headers
  bytes_read=fread(&fbin1,1,sizeof(struct filterbank),infile1);
  bytes_read=fread(&fbin2,1,sizeof(struct filterbank),infile2);

  printf("%f %f\n",fbin1.freq,fbin2.freq);
  printf("%f %f\n",fbin1.bw,fbin2.bw);
  printf("%d %d\n",fbin1.npol,fbin2.npol);
  printf("%d %d\n",fbin1.nchan,fbin2.nchan);
  printf("%f %f\n",fbin1.fsamp,fbin2.fsamp);
  printf("%f %f\n",fbin1.tsamp,fbin2.tsamp);

  // Close
  fclose(infile1);
  fclose(infile2);
  fclose(outfile);

  // Free
  fftwf_free(rp1);
  fftwf_free(rp2);
  fftwf_free(cp1);
  fftwf_free(cp2);

  return 0;
}
