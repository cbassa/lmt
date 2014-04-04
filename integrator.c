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

// This module integrates (sums) complex amplitudes into subints and folds them in phase, leaving the same number of frequency and polarisations channels as its input; it is hardwired to expect 8 polarisations; there is a choice of two if clauses below, which allows or forbids the inclusion of incomplete subints; the number of phase bins, nbin, can be passed in as an option, and the programme will attempt to calculate a suitable number if it is not; if nbin==1 then no folding is performed

int main(int argc,char *argv[])
{
  unsigned int nbin=0,nmax=1024,nsamp=128; // If a value of nbin is passed in, it gets used as long as nbin>=1; if not, default nbin>=1 gets used or default nbin<1 makes the programme calculate a value for nbin based on pulsar period, pulsar DM, tsamp and nmax (the calculated value will be a power of 2 less than or equal to nmax, so nmax may as well be a power of 2); nsamp is number of samples folded into each subint
  int *bintally,arg,warncount,vals_read,binno,sampcount,i,j,nel;
  float *bintally_float;
  // Pulsar phase as fraction of a rotation, pulsar phase as bin number, pulsar period in seconds, sampling interval in bins, pulsar dispersion measure, other stuff
  double phase_start,bin_start,period_start,tsamp_bins,dm,signf,minf,nextf,tmax,lg2=0.69314718055994530941723212145818;
  // FFT in each 'polarisation' (acquired from correlator), and integrated, folded FFTs to output (complex, nchan frequency channels, subints of nsamp samples each, nbin profile bins where nbin==1 for a calibrator)
  fftwf_complex *rp1,*rp2,*rp3,*rp4,*rp5,*rp6,*rp7,*rp8,*ip1,*ip2,*ip3,*ip4,*ip5,*ip6,*ip7,*ip8;
  // File names (these may be FIFOs or actual files)
  FILE *infile,*outfile,*parfile;
  char infname[LIM],outfname[LIM],parfname[LIM];
  struct filterbank fbin,fbout;

  // Decode options: i=input file name (required), o=output file name (required), p=parfile name (required if nbin!=1), t=number of samples per subint (optional; default is preset above), b=number of bins in a folded profile (optional; default is set below according to pulsar period, pulsar DM and tsamp)
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
      nbin=(unsigned int) atoi(optarg); // If a value of nbin is passed in, it overrides any calculation of the value as long as it is greater than zero
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
  vals_read=fread(&fbin,1,sizeof(struct filterbank),infile);
  if (vals_read<1) {
    fprintf(stderr,"Error reading header\n");
    exit;
  }

  // Copy filterbank struct
  fbout.mjd_start=fbin.mjd_start;
  fbout.intmjd=fbin.intmjd;
  fbout.intsec=fbin.intsec;
  strcpy(fbout.source,fbin.source);
  strcpy(fbout.telescope,fbin.telescope);
  strcpy(fbout.instrument,fbin.instrument);
  fbout.freq=fbin.freq;
  fbout.bw=fbin.bw; // This is nchan*fsamp
  fbout.npol=fbin.npol; // Retain all polarisation channels
  fbout.nbit=fbin.nbit; // Number of bits per value (size of float)
  fbout.ndim=fbin.ndim; // Should be complex input
  fbout.nchan=fbin.nchan; // Number of channels
  fbout.fsamp=fbin.fsamp; // Channel size
  fbout.tsamp=fbin.tsamp*nsamp; // Updated sample size

  // Unless no folding is required, we need to find out about the pulsar
  if (nbin!=1) {
    // Try to open par file and fail if it doesn't exist
    parfile=fopen(parfname,"r");
    if (parfile==NULL) {
      fprintf(stderr,"Error opening %s\n",parfname);
      exit;
    }
    // Get pulsar phase and period at initial MJD of dada file using a par file, so we know which phase bins to fold samples into (converted integer part of MJD into an int because it is passed in as an unsigned int; telescope site is currently hardwired to "h" for Effelsberg; source name is passed without initial B or J; parfname is relative path to par file, including directory structure and par file name)
    predict((int)fbin.intmjd,fbin.mjd_start-(double)fbin.intmjd,fbin.source+1,parfname,"h",&phase_start,&period_start,&dm);
    // Use period and DM to calculate a suitable number of bins if it has not been provided (as long as nmax is a sensible number)
    if (nbin<1 && nmax>1) {
      signf=fbin.freq>=0.0?1.0:-1.0; // Just in case centre frequency is negative (which it probably shouldn't be)
      minf=fbin.freq-signf*fabs(fbin.bw)/2.0; // Lowest frequency (closest to zero) in band, assuming band doesn't cross zero (bottom of lowest channel)
      nextf=minf+signf*fabs(fbin.fsamp); // Next frequency away from zero after lowest (top of lowest channel)
      tmax=fabs(4148.80642*dm*(1.0/minf/minf-1.0/nextf/nextf)); // Dispersion smearing in lowest channel if dedispersion is incoherent (assumes times are in units of seconds, frequencies in MHz and DM in cm^-3 pc; the fabs function is just in case it somehows comes out negative, which it shouldn't unless band crosses zero)
      tmax=period_start/(tmax>fbin.tsamp?tmax:fbin.tsamp); // tmax becomes the maximum number of bins using the greater of dispersion smearing and sampling time
      nbin=(int)pow(2.0,floor(log((double)nmax<tmax?(double)nmax:tmax)/lg2)); // Maximum number of bins allowed, using the lesser of tmax and nmax and then going to the next integer power of 2 less than or equal to that
    }
    // Use period and phase for folding if a number of bins greater than 1 is to be used
    if (nbin>1) { // Where a sensible nbin is provided, it doesn't have to be a power of 2
      }
      // Get bin number at initial MJD of dada file
      bin_start=(double)nbin*phase_start;
      // Get sampling interval in units of bins
      tsamp_bins=(double)nbin*fbin.tsamp/period_start;
      printf("Integrator: folding with %d phase bins\n",nbin);
      // Warn if there is sub-sample folding, but allow it if nbin is provided
      if (tsamp_bins>1) {
	printf("Integrator warning: phase bin interval is less than sampling interval\n");
    }
  }
  // Don't use pulsar phase and period if 1 bin is to be used, as folding is not needed (calibrator or pulsar timeseries)
  if (nbin<=1) { // If, after everything, nbin still comes out as a silly number, just make it 1
    printf("Integrator: no folding performed\n");
    bin_start=0;
    tsamp_bins=1;
  }

  // Number of complex elements in each rp and ip array
  nel=fbout.nchan*nbin;

  // Print information
  printf("Integrator: integrating %d spectra into each subint, giving %g us sampling\n",nsamp,fbout.tsamp*1e6);
  printf("Integrator: converting to %d polarizations, %d bits per value\n",fbout.npol,fbout.nbit); // Currently these values don't change from input to output

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

  // Count of total samples processed
  sampcount=0;
  // A loop that keeps going until it's broken (see break statements later)
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
      // Read data
      fread(rp1,sizeof(fftwf_complex),fbin.nchan,infile);
      fread(rp2,sizeof(fftwf_complex),fbin.nchan,infile);
      fread(rp3,sizeof(fftwf_complex),fbin.nchan,infile);
      fread(rp4,sizeof(fftwf_complex),fbin.nchan,infile);
      fread(rp5,sizeof(fftwf_complex),fbin.nchan,infile);
      fread(rp6,sizeof(fftwf_complex),fbin.nchan,infile);
      fread(rp7,sizeof(fftwf_complex),fbin.nchan,infile);
      vals_read=fread(rp8,sizeof(fftwf_complex),fbin.nchan,infile);
      // Exit subint loop when a complete spectrum cannot be read (this shouldn't normally happen until no spectrum at all can be read), so incomplete spectra are not written out and we can decide what to do with incomplete subints later (see next break statement)
      if (vals_read<fbin.nchan) {
	if (vals_read>0)
	  printf("Integrator warning: read incomplete spectrum at end of file, with only %d values in each polarisation instead of %d; these values will be discarded and not integrated.\n",vals_read,fbin.nchan);
	break;

      // Profile bin into which each value must be folded (the first part of the calculation gives a bin number including fractional part, but the floor function rounds this down to an integer, before the modulo operator is applied to give the remainder from division by nbin; so bin 0, for example, is home to everything from 0<=bin<1, rather than, say, -0.5<=bin<0.5)
      binno=floor((double)sampcount*tsamp_bins+bin_start);
      binno%=nbin;
      // Keep track of number of values going into each bin
      bintally[binno]++;

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

    // Exit file loop when a complete spectrum cannot be read, so a partial subint is not written out at the end (first if statement), OR allow partial subints but not empty subints to be written out (second if statement including else-if)
    //    if (vals_read<fbin.nchan) {
    //      if (i>0)
    //	      printf("Integrator warning: discarded %d spectra in each polarisation at the end of the file, because they did not form a complete subint.\n",nsamp-i);
    //      break;
    //    }
    if (i==0) {
      break;
    }
    else if (i<nsamp) {
      printf("Integrator warning: integrated %d spectra in each polarisation at the end of the file, even though they did not form a complete subint\n",i);
    }

    for(i=0;i<=nbin;i++) {
      if (bintally[i]>0) {
	bintally_float[i]=(float)bintally[i];
      }
      else {
	bintally_float[i]=1.0;
	if (warncount==0) {
	  warncount=1;
	  printf("Integrator warning: some profile bins are empty\n");
	}
      }
    }
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
