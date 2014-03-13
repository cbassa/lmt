#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <getopt.h>
#include "filterbank.h"

#define LIM 256

struct filterbank initialize_visibility(struct filterbank fb1,struct filterbank fb2)
{
  double fmin1,fmax1,fmin2,fmax2,fmin,fmax;
  struct filterbank fb;

  // Copy other parameters
  fb.mjd_start=fb1.mjd_start;
  fb.intmjd=fb1.intmjd;
  fb.intsec=fb1.intsec;
  strcpy(fb.source,fb1.source);
  strcpy(fb.telescope,fb1.telescope);
  strcpy(fb.instrument,fb1.instrument);
  fb.tsamp=fb1.tsamp;
  fb.fsamp=fb1.fsamp;
  fb.nbit=fb1.nbit;
  fb.ndim=fb1.ndim;
  fb.npol=fb1.npol;

  // Determine frequency overlap
  fmin1=fb1.freq-0.5*fabs(fb1.bw);
  fmax1=fb1.freq+0.5*fabs(fb1.bw);
  fmin2=fb2.freq-0.5*fabs(fb2.bw);
  fmax2=fb2.freq+0.5*fabs(fb2.bw);

  if (fmin1>=fmin2)
    fmin=fmin1;
  else
    fmin=fmin2;
  if (fmax1<=fmax2)
    fmax=fmax1;
  else
    fmax=fmax2;

  // Center frequency
  fb.freq=0.5*(fmax+fmin);
  fb.bw=fmax-fmin;

  // Number of channels
  fb.nchan=(int) ceil(fb.bw/fb.fsamp);

  printf("Correlator: center frequency: %g MHz, bandwidth: %g MHz, %d channels\n",fb.freq,fb.bw,fb.nchan);

  return fb;
}

int main(int argc,char *argv[])
{
  int arg=0;
  int i=0,i1,i2,j,k=0;
  fftwf_complex *c1p1,*c1p2,*c2p1,*c2p2,*ccp1,*ccp2,*cdp1,*cdp2,*a1p1,*a1p2,*a2p1,*a2p2;
  FILE *infile1,*infile2,*outfile;
  char infname1[LIM],infname2[LIM],outfname[LIM];
  int bytes_read;
  struct filterbank fbin1,fbin2,fbout;
  double freq,freq1,freq2;

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

  // Generate output struct, compute frequency overlap, etc
  fbout=initialize_visibility(fbin1,fbin2);

    // Write header struct
  fwrite(&fbout,1,sizeof(struct filterbank),outfile);

  // Allocate buffers
  c1p1=fftwf_malloc(sizeof(fftwf_complex)*fbin1.nchan);
  c1p2=fftwf_malloc(sizeof(fftwf_complex)*fbin1.nchan);
  c2p1=fftwf_malloc(sizeof(fftwf_complex)*fbin2.nchan);
  c2p2=fftwf_malloc(sizeof(fftwf_complex)*fbin2.nchan);
  ccp1=fftwf_malloc(sizeof(fftwf_complex)*fbout.nchan);
  ccp2=fftwf_malloc(sizeof(fftwf_complex)*fbout.nchan);
  cdp1=fftwf_malloc(sizeof(fftwf_complex)*fbout.nchan);
  cdp2=fftwf_malloc(sizeof(fftwf_complex)*fbout.nchan);
  a1p1=fftwf_malloc(sizeof(fftwf_complex)*fbout.nchan);
  a1p2=fftwf_malloc(sizeof(fftwf_complex)*fbout.nchan);
  a2p1=fftwf_malloc(sizeof(fftwf_complex)*fbout.nchan);
  a2p2=fftwf_malloc(sizeof(fftwf_complex)*fbout.nchan);

  for (;;) {
    // Read buffers
    bytes_read=fread(c1p1,sizeof(fftwf_complex),fbin1.nchan,infile1);
    bytes_read=fread(c1p2,sizeof(fftwf_complex),fbin1.nchan,infile1);

    // Exit when buffer is empty
    if (bytes_read==0)
      break;

    // Read buffers
    bytes_read=fread(c2p1,sizeof(fftwf_complex),fbin2.nchan,infile2);
    bytes_read=fread(c2p2,sizeof(fftwf_complex),fbin2.nchan,infile2);

    // Exit when buffer is empty
    if (bytes_read==0)
      break;


    for (i=0;i<fbout.nchan;i++) {
      // Compute frequencies
      freq=fbout.freq-0.5*fbout.bw+fbout.bw*(float) (i+0.5)/(float) fbout.nchan;

      // Channel numbers
      i1=(int) fbin1.nchan*(freq-fbin1.freq+0.5*fbin1.bw)/fbin1.bw;
      i2=(int) fbin2.nchan*(freq-fbin2.freq+0.5*fbin2.bw)/fbin2.bw;
      
      //      r=r1*r2+i1*i2;
      //      i=i1*r2-r1*i2;
      // Auto correlation 1
      a1p1[i][0]=c1p1[i][0]*c1p1[i][0]+c1p1[i][1]*c1p1[i][1];
      a1p1[i][1]=c1p1[i][1]*c1p1[i][0]-c1p1[i][0]*c1p1[i][1];
      a1p2[i][0]=c1p2[i][0]*c1p2[i][0]+c1p2[i][1]*c1p2[i][1];
      a1p2[i][1]=c1p2[i][1]*c1p2[i][0]-c1p2[i][0]*c1p2[i][1];

      // Auto correlation 2
      a2p1[i][0]=c2p1[i][0]*c2p1[i][0]+c2p1[i][1]*c2p1[i][1];
      a2p1[i][1]=c2p1[i][1]*c2p1[i][0]-c2p1[i][0]*c2p1[i][1];
      a2p2[i][0]=c2p2[i][0]*c2p2[i][0]+c2p2[i][1]*c2p2[i][1];
      a2p2[i][1]=c2p2[i][1]*c2p2[i][0]-c2p2[i][0]*c2p2[i][1];

      // Cross correlation L1L2 R1R2
      cdp1[i][0]=c1p1[i][0]*c2p1[i][0]+c1p1[i][1]*c2p1[i][1];
      cdp1[i][1]=c1p1[i][1]*c2p1[i][0]-c1p1[i][0]*c2p1[i][1];
      cdp2[i][0]=c1p2[i][0]*c2p2[i][0]+c1p2[i][1]*c2p2[i][1];
      cdp2[i][1]=c1p2[i][1]*c2p2[i][0]-c1p2[i][0]*c2p2[i][1];
      // Cross correlation L1R2 R1L2
      ccp1[i][0]=c1p1[i][0]*c2p2[i][0]+c1p1[i][1]*c2p1[i][1];
      ccp1[i][1]=c1p1[i][1]*c2p2[i][0]-c1p1[i][0]*c2p1[i][1];
      ccp2[i][0]=c1p2[i][0]*c2p1[i][0]+c1p2[i][1]*c1p2[i][1];
      ccp2[i][1]=c1p2[i][1]*c2p1[i][0]-c1p2[i][0]*c1p2[i][1];
    }
  }

  // Close
  fclose(infile1);
  fclose(infile2);
  fclose(outfile);

  // Free
  fftwf_free(c1p1);
  fftwf_free(c1p2);
  fftwf_free(c2p1);
  fftwf_free(c2p2);
  fftwf_free(ccp1);
  fftwf_free(ccp2);
  fftwf_free(cdp1);
  fftwf_free(cdp2);
  fftwf_free(a1p1);
  fftwf_free(a1p2);
  fftwf_free(a2p1);
  fftwf_free(a2p2);

  return 0;
}
