#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <getopt.h>
#include "filterbank.h"
#include "lib/delays.h"

#define LIM 256

void Usage()
{
  printf("Usage: reader -i <input file> -o <output file> [-b <blocksize> (64000)] -c <calcfile> -f <\"fracdelay_pol1 fracdelay_pol2\"> -p <\"phaseoffset_pol1 phaseoffset_pol2\">\n");
}

int main(int argc,char *argv[])
{
  int arg=0;
  int i=0,j,k=0;
  fftwf_complex *rp1,*rp2,*cp1,*cp2;
  FILE *infile,*outfile;
  char infname[LIM],outfname[LIM];
  int bytes_read;
  calc_type calc;
  phases_type phases;
  struct filterbank fbin,fbout;
  int iset=0, oset=0, cset=0, fset=0, pset=0;

  // Decode options
  while ((arg=getopt(argc,argv,"i:o:b:c:f:p:"))!=-1) {
    switch (arg) {

    case 'i':
      strcpy(infname,optarg);
      iset=1;
      break;

    case 'o':
      strcpy(outfname,optarg);
      oset=1;
      break;

      // Get the blocksize
    case 'b':
      blocksize=(unsigned int) atoi(optarg);
      break;

      // Get calc filename 
    case 'c':
      calc.filename = optarg;
      cset=1;
      break;

      // Get fractional delays for pol1 and pol2 
    case 'f':
      sprintf(optarg, "%f %f\n", &phases.frac_delay_pol1, &phases.frac_delay_pol1);
      fset=1;
      break;

      // Get phase offsets for pol1 and pol2 
    case 'p':
      sprintf(optarg, "%f %f\n", &phases.phase_offset_pol1, &phases.phase_offset_pol2);
      pset=1;
      break;

    default:
      return 0;
    }
  }

  // Check that all necessary parameters are set.
  if (!iset){ fprintf(stderr, "Please provide input parameter with -i\n"); Usage(); exit(0);}
  if (!oset){ fprintf(stderr, "Please provide output parameter with -o\n"); Usage(); exit(0);}
  if (!cset){ fprintf(stderr, "Please provide calcfilename with -c\n"); Usage(); exit(0);}
  if (!fset){ fprintf(stderr, "Please provide fraction delays with -f\n"); Usage(); exit(0);}
  if (!pset){ fprintf(stderr, "Please provide phase offsets with -p\n"); Usage(); exit(0);}

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
  fbout.npol=fbin.npol;
  fbout.nbit=fbin.nbit;
  fbout.ndim=fbin.ndim;
  fbout.nchan=fbin.nchan;
  fbout.fsamp=fbin.fsamp;
  fbout.tsamp=fbin.tsamp;

  // Print information
  printf("Donothing_fb: just copying a filterbank struct\n");

  // Write header struct
  fwrite(&fbout,1,sizeof(struct filterbank),outfile);

  // Allocate buffers
  rp1=fftwf_malloc(sizeof(fftwf_complex)*fbin.nchan);
  rp2=fftwf_malloc(sizeof(fftwf_complex)*fbin.nchan);
  cp1=fftwf_malloc(sizeof(fftwf_complex)*fbout.nchan);
  cp2=fftwf_malloc(sizeof(fftwf_complex)*fbout.nchan);

  // Loop over file contents
  for (;;) {
    // Read buffers
    bytes_read=fread(rp1,sizeof(fftwf_complex),fbin.nchan,infile);
    bytes_read=fread(rp2,sizeof(fftwf_complex),fbin.nchan,infile);

    // Exit when buffer is empty
    if (bytes_read==0)
      break;

    // Copy
    for (j=0;j<fbin.nchan;j++) {
      cp1[j][0]=rp1[j][0];
      cp1[j][1]=rp1[j][1];
      cp2[j][0]=rp2[j][0];
      cp2[j][1]=rp2[j][1];
    }

    // Write
    fwrite(cp1,sizeof(fftwf_complex),fbout.nchan,outfile);
    fwrite(cp2,sizeof(fftwf_complex),fbout.nchan,outfile);
  } 

  // Close
  fclose(infile);
  fclose(outfile);

  // Free
  fftwf_free(rp1);
  fftwf_free(rp2);
  fftwf_free(cp1);
  fftwf_free(cp2);

  return 0;
}
