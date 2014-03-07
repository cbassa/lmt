#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <getopt.h>
#include "filterbank.h"
#include "predict.h"

// File name character limit
#define LIM 256

// This module integrates. predict routine gets period. nbins is passed in. addtovis adds in lcorr. Calibrator option?

int main(int argc,char *argv[])
{
  int arg=0;
  int *bintally,binno,sampcount,i,j,nel,nsamp=128,nbin=128;
  float *bintally_float;
  // Pulsar phase as fraction of a rotation, pulsar phase as bin number, pulsar period in seconds and sampling interval in bins
  double phase_start,bin_start,period_start,tsamp_bins;
  // FFT in each 'polarisation' (acquired from correlator), and integrated, folded FFTs to output (complex, nchan frequency channels, subints of nsamp samples each, nbin profile bins where nbin==1 for a calibrator)
  fftwf_complex *rp1,*rp2,*rp3,*rp4,*rp5,*rp6,*rp7,*rp8,*ip1,*ip2,*ip3,*ip4,*ip5,*ip6,*ip7,*ip8;
  // File names (these may be FIFOs or actual files)
  FILE *infile,*outfile,*parfile;
  char infname[LIM],outfname[LIM],parfname[LIM];
  int bytes_read;
  struct filterbank fbin,fbout;

  // Decode options: i=input file name, o=output file name, p=parfile name, t=number of samples per subint (default is preset above), b=number of bins in a folded profile (default is preset above)
  while ((arg=getopt(argc,argv,"i:o:p:t:b:"))!=-1) {
    switch (arg) {

    case 'i':
      strcpy(infname,optarg);
      break;

    case 'o':
      strcpy(outfname,optarg);
      break;

    case 'p':
      strcpy(parfname,optarg);
      break;

    case 't':
      nsamp=(unsigned int) atoi(optarg);
      break;
      
    case 'b':
      nbin=(unsigned int) atoi(optarg);
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

  // Open par file
  parfile=fopen(parfname,"r");

  // Check if par file exists
  if (parfile==NULL) {
    fprintf(stderr,"Error opening %s\n",parfname);
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
  fbout.bw=fbin.bw; // This should be nchan*fsamp, presumably...?
  fbout.npol=fbin.npol; // Retain all polarisation channels
  fbout.nbit=fbin.nbit; // Number of bits per value (size of float)
  fbout.ndim=fbin.ndim; // Should be complex input
  fbout.nchan=fbin.nchan; // Number of channels
  fbout.fsamp=fbin.fsamp; // Channel size
  fbout.tsamp=fbin.tsamp*nsamp; // Updated sample size

  // Number of complex elements in each rp and ip array
  nel=fbout.nchan*nbin;

  // Print information
  printf("Integrator: integrating %d spectra, giving %g us sampling\n",nsamp,fbout.tsamp*1e6);
  printf("Integrator: converting to %d polarizations, %d bit\n",fbout.npol,fbout.nbit); // Currently these values don't change from input to output

  // Write header struct
  fwrite(&fbout,1,sizeof(struct filterbank),outfile);

  // Allocate buffers
  bintally=malloc(sizeof(int)*nbin);
  bintally_float=malloc(sizeof(float)*nbin);
  rp1=fftwf_malloc(sizeof(fftwf_complex)*fbin.nchan);
  rp2=fftwf_malloc(sizeof(fftwf_complex)*fbin.nchan);
  rp3=fftwf_malloc(sizeof(fftwf_complex)*fbin.nchan);
  rp4=fftwf_malloc(sizeof(fftwf_complex)*fbin.nchan);
  rp5=fftwf_malloc(sizeof(fftwf_complex)*fbin.nchan);
  rp6=fftwf_malloc(sizeof(fftwf_complex)*fbin.nchan);
  rp7=fftwf_malloc(sizeof(fftwf_complex)*fbin.nchan);
  rp8=fftwf_malloc(sizeof(fftwf_complex)*fbin.nchan);
  ip1=fftwf_malloc(sizeof(fftwf_complex)*nel);
  ip2=fftwf_malloc(sizeof(fftwf_complex)*nel);
  ip3=fftwf_malloc(sizeof(fftwf_complex)*nel);
  ip4=fftwf_malloc(sizeof(fftwf_complex)*nel);
  ip5=fftwf_malloc(sizeof(fftwf_complex)*nel);
  ip6=fftwf_malloc(sizeof(fftwf_complex)*nel);
  ip7=fftwf_malloc(sizeof(fftwf_complex)*nel);
  ip8=fftwf_malloc(sizeof(fftwf_complex)*nel);

  if (nbin>1) {
    // Get pulsar phase and period at initial MJD of dada file using a par file, so we know which phase bins to fold samples into (converted integer part of MJD into an int because it is passed in as an unsigned int; telescope site is currently hardwired to "h" for Effelsberg; source name is passed without initial B or J; parfname is relative path to par file, including directory structure and par file name)
    predict((int)fbin.intmjd,fbin.mjd_start-(double)fbin.intmjd,fbin.source+1,parfname,"h",&phase_start,&period_start);
    // Get bin number at initial MJD of dada file
    bin_start=(double)nbin*phase_start;
    // Get sampling interval in units of bins
    tsamp_bins=(double)nbin*fbin.tsamp/period_start;
  }
  // Don't look for pulsar phase and period if nbin=1 (generally for calibrators)
  else {
    bin_start=0;
    tsamp_bins=1;
  }

  // Count of total samples processed
  sampcount=0;
  // Read buffers up to the end
  for (;;) {
    // Reset buffers
    for (j=0;j<nbin;j++)
      bintally[j]=0;
    for (j=0;j<nel;j++) {
      ip1[j][0]=0.0;
      ip1[j][1]=0.0;
      ip2[j][0]=0.0;
      ip2[j][1]=0.0;
      ip3[j][0]=0.0;
      ip3[j][1]=0.0;
      ip4[j][0]=0.0;
      ip4[j][1]=0.0;
      ip5[j][0]=0.0;
      ip5[j][1]=0.0;
      ip6[j][0]=0.0;
      ip6[j][1]=0.0;
      ip7[j][0]=0.0;
      ip7[j][1]=0.0;
      ip8[j][0]=0.0;
      ip8[j][1]=0.0;
    }

    // Accumulate for each subint (and fold if necessary) by reading FFT data from input file (or from another module if using FIFOs)
    for (i=0;i<nsamp;i++) {
      // Profile bin into which each value must be folded (the first part of the calculation gives a bin number including fractional part, but the floor function rounds this down to an integer, before the modulo operator is applied to give the remainder from division by nbin; so bin 0, for example, is home to everything from 0<=bin<1, rather than, say, -0.5<=bin<0.5)
      binno=floor((double)sampcount*tsamp_bins+bin_start);
      binno%=nbin;
      // Keep track of number of values going into each bin
      bintally[binno]++;
      // Read data
      bytes_read=fread(rp1,sizeof(fftwf_complex),fbin.nchan,infile);
      bytes_read=fread(rp2,sizeof(fftwf_complex),fbin.nchan,infile);
      bytes_read=fread(rp3,sizeof(fftwf_complex),fbin.nchan,infile);
      bytes_read=fread(rp4,sizeof(fftwf_complex),fbin.nchan,infile);
      bytes_read=fread(rp5,sizeof(fftwf_complex),fbin.nchan,infile);
      bytes_read=fread(rp6,sizeof(fftwf_complex),fbin.nchan,infile);
      bytes_read=fread(rp7,sizeof(fftwf_complex),fbin.nchan,infile);
      bytes_read=fread(rp8,sizeof(fftwf_complex),fbin.nchan,infile);

      // Fold subint for each channel and polarisation
      for (j=0;j<fbin.nchan;j++) {
	// Bins becomes the first dimension in packing, before channels, polarisations and time
	ip1[binno][0]+=rp1[j][0];
	ip1[binno][1]+=rp1[j][1];
	ip2[binno][0]+=rp2[j][0];
	ip2[binno][1]+=rp2[j][1];
	ip3[binno][0]+=rp3[j][0];
	ip3[binno][1]+=rp3[j][1];
	ip4[binno][0]+=rp4[j][0];
	ip4[binno][1]+=rp4[j][1];
	ip5[binno][0]+=rp5[j][0];
	ip5[binno][1]+=rp5[j][1];
	ip6[binno][0]+=rp6[j][0];
	ip6[binno][1]+=rp6[j][1];
	ip7[binno][0]+=rp7[j][0];
	ip7[binno][1]+=rp7[j][1];
	ip8[binno][0]+=rp8[j][0];
	ip8[binno][1]+=rp8[j][1];
	binno+=nbin;
      }
      // Move sample number along
      sampcount++;
    }

    // Exit when buffer is empty (WARNING: this could go off the end of the file with an incomplete subint)
    if (bytes_read==0)
      break;

    for(i=0;i<=nbin;i++)
      bintally_float[i]=(float)bintally[i];
    // Within each bin, scale bins by a factor of the number of samples in that bin
    for (j=0,binno=0;j<fbin.nchan;j++) {
      for (i=0;i<nbin;i++) {
	ip1[binno][0]/=bintally_float[i];
	ip1[binno][1]/=bintally_float[i];
	ip2[binno][0]/=bintally_float[i];
	ip2[binno][1]/=bintally_float[i];
	ip3[binno][0]/=bintally_float[i];
	ip3[binno][1]/=bintally_float[i];
	ip4[binno][0]/=bintally_float[i];
	ip4[binno][1]/=bintally_float[i];
	ip5[binno][0]/=bintally_float[i];
	ip5[binno][1]/=bintally_float[i];
	ip6[binno][0]/=bintally_float[i];
	ip6[binno][1]/=bintally_float[i];
	ip7[binno][0]/=bintally_float[i];
	ip7[binno][1]/=bintally_float[i];
	ip8[binno][0]/=bintally_float[i];
	ip8[binno][1]/=bintally_float[i];
	binno++;
      }
    }

    // Write out the integrated power in each channel for each subint
    fwrite(ip1,sizeof(fftw_complex),fbin.nchan,outfile);
    fwrite(ip2,sizeof(fftw_complex),fbin.nchan,outfile);
    fwrite(ip3,sizeof(fftw_complex),fbin.nchan,outfile);
    fwrite(ip4,sizeof(fftw_complex),fbin.nchan,outfile);
    fwrite(ip5,sizeof(fftw_complex),fbin.nchan,outfile);
    fwrite(ip6,sizeof(fftw_complex),fbin.nchan,outfile);
    fwrite(ip7,sizeof(fftw_complex),fbin.nchan,outfile);
    fwrite(ip8,sizeof(fftw_complex),fbin.nchan,outfile);
  } 

  // Close
  fclose(infile);
  fclose(outfile);

  // Free
  fftwf_free(rp1);
  fftwf_free(rp2);
  fftwf_free(rp3);
  fftwf_free(rp4);
  fftwf_free(rp5);
  fftwf_free(rp6);
  fftwf_free(rp7);
  fftwf_free(rp8);
  free(ip1);
  free(ip2);
  free(ip3);
  free(ip4);
  free(ip5);
  free(ip6);
  free(ip7);
  free(ip8);
  free(bintally);
  free(bintally_float);

  return 0;
}
