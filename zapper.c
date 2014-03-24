#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <getopt.h>
#include "filterbank.h"
#include "gsl/gsl_sf_erf.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_integration.h"

// File name character limit
#define LIM 256

// This module uses spectral kurtosis to create a mask of good (1) and bad (0) values corresponding to the channelised data. It analyses total power and therefore masks both polarisations (or neither) at each time and frequency. The temporal and spectral masking resolutions are given by nskz*nav and nchan (you can average nav spectra together if you want), but this does not affect the resolutions of the data. It does not mask any data in channel 0 (or channel nbin, if present after the FFT of a real time series).

int sk_thresh5(int N_int, int M_int, float s_float, float d_float, float sk_lims_float[]);

int main(int argc,char *argv[])
{
  // nskz is called M in the Nita & Gary papers (number of blocks of nav FFTs of each polarisation used to compute each SK value), while nav is called N (number of spectra averaged together before calculating each SK value)
  int nskz=2048,nav=1; // e.g. 2048 and 1 (default values)
  // Options: 1/0 for per-channel detection or not; 1/0 to normalise summed power in each channel and each polarisation by total bandwidth power in that polarisation, or not (only applies to noscr option; using it may make you miss RFI, but is necessary if you have a very strong signal which has smooth time variation across nskz*nav*nchan time samples and which you want to let through); 1/0 for frequency-scrunched detection (across full bandwidth) or not
  int noscr=1,norm_noscr=0,fscr=0; // e.g. 1, 0 and 0 (default values)
  // Probabilities of falsely excluding SK values from non-RFI signals in channelised and frequency-scrunched data, in terms of sigma; this is used to set SK thresholds, e.g. if sig_noscr=3.0 you will falsely exclude about 0.27% of good data through channelised zapping; if you set sig_noscr or sig_fscr low you will exclude more good data, but if you set them high you may include more RFI
  float sig_noscr=3.0,sig_fscr=3.0; // e.g. 3.0 and 3.0 (default values)
  // Desired spectral resolution of masking
  double sk_fsamp=0.0; // Default value of 0 makes programme assign sk_fsamp equal to fsamp (existing frequency resolution) unless another value is passed in
  // Shape of expected distribution of most SK values (e.g. 1.0 for standard exponential distribution, corresponding to power made from a complex value such as an FFT amplitude)
  float d=1.0;

  int *k,count_chunk,count_noscr,count_fscr,i,j,l,iav,nchan,ichan,ndim,bytes_read,nsamp,nav2,nskz_fscr,chan_fac,sk_lim_status;
  char *mask,infname[LIM],outfname[LIM];
  fftwf_complex *rp1,*rp2; // Filterbank data
  FILE *infile,*outfile;
  float *pp1,*s1p1,*s2p1,*skp1,sk_lims[2],sk_lims_fscr[2],nchan_float,ppnosum,ppsum1_noscr,apn1,s1p1_fscr,s2p1_fscr,skp1_fscr,zap_noscr,zap_fscr,Nd;
  struct filterbank fbin,fbout;

  // Assign some useful values
  nav2=nav*2; // nav is doubled because we are using total power made from two polarisations
  ppnosum=(float)nav2*(float)nskz;

  // Decode options: i=input file name (required), o=output file name (required), b=desired masking frequency resolution (optional; default is fbin.fsamp), m=
  while ((arg=getopt(argc,argv,"i:o:b:m:n:c:s:r:f:g:"))!=-1) {
    switch (arg) {

    case 'i':
      strcpy(infname,optarg);
      break;

    case 'o':
      strcpy(outfname,optarg);
      break;

    case 'b':
      sk_fsamp=strtod(optarg,NULL); // Desired frequency resolution in MHz (called 'b' for bandwidth) of RFI masking (relevant if -c 1 is used)
      break;

    case 'm':
      nskz=(unsigned int) atoi(optarg);
      break;

    case 'n':
      nav=(unsigned int) atoi(optarg);
      break;

    case 'c':
      noscr=(unsigned int) atoi(optarg); // 1 or 0
      break;

    case 's':
      sig_noscr=strtof(optarg,NULL);
      break;

    case 'r':
      norm_noscr=(unsigned int) atoi(optarg); // 1 or 0 (1 to normalise or rescale - hence the 'r' - values by the total bandwidth power if -c 1 is used)
      break;

    case 'f':
      fscr=(unsigned int) atoi(optarg); // 1 or 0
      break;

    case 'g':
      sig_fscr=(float) strtof(optarg,NULL);
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
  fbout.bw=fbin.bw; // This is nchan*fsamp
  fbout.npol=fbin.npol; // Retain all polarisation channels
  fbout.nbit=fbin.nbit; // Number of bits per value (size of float)
  fbout.ndim=fbin.ndim; // Should be complex input
  fbout.nchan=fbin.nchan; // Number of channels
  fbout.fsamp=fbin.fsamp; // Channel size
  fbout.tsamp=fbin.tsamp; // Sample size
  // Write header struct
  fwrite(&fbout,1,sizeof(struct filterbank),outfile);
  
  // Allocate memory
  rp1=fftwf_malloc(sizeof(fftwf_complex)*fbin.nchan);
  rp2=fftwf_malloc(sizeof(fftwf_complex)*fbin.nchan);
  mask=fftwf_malloc(sizeof(char)*nchan);

  // Choose an integer factor by which to frequency-scrunch when masking, which is the integer nearest to the desired value according to sk_fsamp
  chan_fac=(int)round(sk_fsamp/fbin.fsamp);
  if (chan_fac<1)
    chan_fac=1; // Frequency-scrunching factor must be at least 1
  if (chan_fac>nchan)
    chan_fac=nchan; // Frequency-scrunching factor must be at most nchan
  nchan=(int)ceil((double)fbin.nchan/(double)chan_fac); // Number of frequency channels for masking, including a part-channel at the end of the band if necessary
  sk_fsamp=fbin.fsamp*(double)chan_fac; // Frequency resolution (channel bandwidth) of masking (except for any part-channel at the end of the band, which will have a smaller value)
  nchan_float=(float)nchan;

  if (noscr==1)
  {
    count_noscr=0;
    s1p1=(float *) malloc(sizeof(float)*nchan);
    s2p1=(float *) malloc(sizeof(float)*nchan);
    skp1=(float *) malloc(sizeof(float)*nchan);
    // Generate SK thresholds using sk_thresh5 function
    sk_lim_status=sk_thresh5(nav2,nskz,sig_noscr,d,sk_lims);
    printf("Zapper: non-scrunched thresholds: %f %f\n",sk_lims[0],sk_lims[1]);
  }
  else
  {
    norm_noscr=0;
  }
  if (fscr==1)
  {
    count_fscr=0;
    // Note that frequency-scrunched values come from the statistics of the set of power values in each channel (except channel 0) and at each time, rather than from a set of power values summed over the whole band at each time (hence the factor of nchan-1 below); the same is done if the requested frequency resolution of masking requires channels to be scrunched together - this is different to how the two polarisation channels are treated, where the power for each value is made from both polarisations summed together, rather than using a set containing values from each polarisation separately; the approach with frequency is to avoid losing sensitivity by averaging too many values together, while the approach with power is to diluting polarised RFI and making it more difficult to detect
    nskz_fscr=nskz*(nchan-1);
    // Generate f-scrunched SK thresholds
    sk_lim_status=sk_thresh5(nav2,nskz_fscr,sig_fscr,d,sk_lims_fscr);
    printf("Zapper: frequency-scrunched thresholds: %f %f\n",sk_lims_fscr[0],sk_lims_fscr[1]);
  }
  Nd=nav2*d;

  // Read first chunk
  bytes_read=fread(buffer,1,buffersize,file);

  count_chunk=0;
  // RFI-zap current chunk and read the next
  while (bytes_read==buffersize) // Don't read/write the last part-chunk of the dada file (so filterbank file may be shorter in time)
  {
    if (bytes_read<buffersize)
      break;

    if (noscr==1)
    {
      for (ichan=0;ichan<nchan;ichan++)
      {
	s1p1[ichan]=0.0;
	s2p1[ichan]=0.0;
      }
    }
    if (fscr==1)
    {
      s1p1_fscr=0.0;
      s2p1_fscr=0.0;
    }

    for (l=0;l<buffersize_power;l++)
      buffer_power[l]=0.0;

    // Loop over nskz blocks of nav summed FFTs in each polarisation
    for (i=0,j=0,l=0;i<buffersize_power;i+=nchan)
    {
      // Sum nav FFTs
      for (iav=0;iav<nav;iav++)
      {  
	// Unpack time series buffer
	if (ndim==2)
	{
	  for (ichan=0;ichan<nchan;ichan++)
	  {
	    rp1[ichan][0]=(float)buffer[j];
	    rp1[ichan][1]=(float)buffer[j+1];
	    rp2[ichan][0]=(float)buffer[j+2];
	    rp2[ichan][1]=(float)buffer[j+3];
	    j+=4;
	  }
	}
	else
	{
	  for (ichan=0;ichan<nsamp;ichan++)
	  {
	    rrp1[ichan]=(float)buffer[j];
	    rrp2[ichan]=(float)buffer[j+1];
	    j+=2;
	  }
	}
	// Forward FFT
	fftwf_execute(ftp1);
	fftwf_execute(ftp2);
	
	// Filterbank power calculated from FFTs (channel 0 is not included in the case of real time series)
	for (ichan=1;ichan<nchan;ichan++)
	  buffer_power[i+k[ichan]]+=cp1[ichan][0]*cp1[ichan][0]+cp1[ichan][1]*cp1[ichan][1]+cp2[ichan][0]*cp2[ichan][0]+cp2[ichan][1]*cp2[ichan][1];
      }
	
      // Everything after this point deals with summed spectra
      if (noscr==1)
      {
	// Normalise or don't normalise power in each channel of summed spectrum by total bandwidth power
	if (norm_noscr==1)
	{
	  if (nchan!=1)
	  {
	    ppsum1_noscr=0.0;
	    // Note that channel 0 is excluded from total bandwidth power
	    for (ichan=1;ichan<nchan;ichan++)
	      ppsum1_noscr+=buffer_power[i+k[ichan]];
	  }
	  else
	  {
	    ppsum1_noscr=buffer_power[i];
	  }
	}
        else
        {
	  ppsum1_noscr=ppnosum;
	}
	  
	// Compute values for spectral kurtosis
	// Channel 0 not used here, as it is not generally comparable between real and complex time series
	for (ichan=1;ichan<nchan;ichan++)
	{
	  apn1=buffer_power[i+k[ichan]]/ppsum1_noscr;
	  s1p1[ichan]+=apn1;
	  s2p1[ichan]+=apn1*apn1;
	}
      }
      if (fscr==1)
      {
	// Channel 0 not used here as it has a different distribution for real time series
	for (ichan=1;ichan<nchan;ichan++)
	{
	  s1p1_fscr+=buffer_power[i+k[ichan]]/nchan_float/ppnosum;
	  s2p1_fscr+=buffer_power[i+k[ichan]]*buffer_power[i+k[ichan]]/nchan_float/nchan_float/ppnosum/ppnosum;
	}
      }
    }

    if (noscr==1)
    {
      // Compute SK sums and make RFI masks (repacked, all values kept for one plot)
      // RFI detection not implemented in channel 0 as it is not generally comparable between real and complex time series; also not implemented in channel nchan, which exists for real time series only
      for (ichan=1;ichan<nchan;ichan++)
      {
	skp1[ichan]=(nskz*Nd+1.0)/(nskz-1.0)*(nskz*s2p1[ichan]/s1p1[ichan]/s1p1[ichan]-1.0);
	// RFI mask for values which exceed SK thresholds (repacked, all values kept for one plot)
	if (skp1[ichan]<sk_lims[0] || skp1[ichan]>sk_lims[1])
	{
	  for (i=k[ichan];i<buffersize_power;i+=nchan)
	    buffer_power[i]=0.0; // Power set to zero in bad channels
	  count_noscr++;
	}
      }
    }
    if (fscr==1)
    {
      skp1_fscr=(nskz_fscr*Nd+1.0)/(nskz_fscr-1.0)*(nskz_fscr*s2p1_fscr/s1p1_fscr/s1p1_fscr-1.0);
      // All channels except 0 and (for real time series) nchan used for frequency-scrunched RFI detection
      if (skp1_fscr<sk_lims_fscr[0] || skp1_fscr>sk_lims_fscr[1])
      {
	for (i=0;i<buffersize_power;i+=nchan)
	{
	  for (ichan=1;ichan<nchan;ichan++)
	    buffer_power[i+k[ichan]]=0.0;
	}
	count_fscr++;
      }
    }

    // Write chunk of filterbank file
    fwrite(buffer_power,sizeof(float),buffersize_power,file_new);
    count_chunk++;

    // Read next chunk
    bytes_read=fread(buffer,1,buffersize,file);
  }

  // Ignore last part-chunk if it exists, because we can't RFI-zap it
  if (bytes_read>0)
    sprintf(end_str," and final %d bytes of time series ",bytes_read);

  // Report zapped percentage
  if (noscr==1)
  {
    zap_noscr=(float)count_noscr/(nchan_float-1.0)*100.0/(float)count_chunk;
    printf("Zapper: %.2f%% of data replaced by channel-zapping (%s%sexcluded from statistics)\n",zap_noscr,chan_str,end_str);
  }
  if (fscr==1)
  {
    zap_fscr=(float)count_fscr/(float)count_chunk*100.0;
    printf("Zapper: %.2f%% of data replaced by bandwidth-zapping (%s%sexcluded from statistics)\n",zap_fscr,chan_str,end_str);
    if (noscr==1)
      printf("Data zapped by both processes are included in both percentages\n");
  }

  // Close files
  fclose(file);
  fclose(file_new);

  // Free
  free(k);
  free(header);
  free(buffer);
  free(buffer_power);
  if (ndim==2)
  {
    fftwf_free(rp1);  
    fftwf_free(rp2);
  }
  else
  {
    free(rrp1);
    free(rrp2);
  }
  fftwf_free(cp1);
  fftwf_free(cp2);
  fftwf_destroy_plan(ftp1);
  fftwf_destroy_plan(ftp2);
  if (noscr==1)
  {
    free(s1p1);
    free(s2p1);
    free(skp1);
  }

  return 0;
}

int sk_thresh5(int N_int, int M_int, float s_float, float d_float, float sk_lims_float[])
{
  /* This function takes two integers, two floats and a 2-element float array as input, fills the array and gives an integer as output */
  /* The function follows PASP papers Nita, Gary & Liu (2007) and Nita & Gary (2010) and MNRAS letter Nita & Gary (2010). It calculates upper and lower thresholds for the spectral kurtosis (SK) RFI estimator. N is the number of FFTs being averaged into power spectra (minimum value 1, corresponding to Nyquist-sampled data with no averaging) (N is NOT the number of samples per FFT); M is number of blocks of N FFTs being used to produce each value of the SK estimator (minimum value 3, otherwise equations don't work); s is the sigma limit of the thresholds; d gives the type of distribution expected from non-RFI (1.0 for standard exponential distribution); sk_lims is the array into which the calculated threshold values will be placed. */
  
  if (M_int<2 || N_int<1)
  {
    printf("Thresholds can't be calculated for your chosen parameters! M must be at least 2 and N must be at least 1. Thresholds have been returned as zero.\n");
    sk_lims_float[0] = 0.0;
    sk_lims_float[1] = 0.0;
    return 1;
  }
  
  else
  {
    /* Initial hard-wired variables, which you might want to change but which seem to work well */
    size_t subdiv_lim = 100; /* Maximum number of intervals into which integrations will be subdivided (this can be quite low for a smoothly varying function such as the probability density function we are integrating, but may need to be larger for large M or N) */
    double Ptol = 0.00001; /* Permitted absolute error in cumulative probability values calculated by integration (this can be small for our smoothly varying function) */
    int max_its = 100; /* Maximum number of numerical integrations allowed when trying to find an initial or best threshold value */
    
    /* Other variables declared */
    int status, ul, fill_lo, fill_hi, n;
    double pi = 3.1415926535897932384626433832795, s=(double)s_float, d=(double)d_float, Nd = (double)N_int*d, M = (double)M_int, sk_lims[2], NN, NN1, M1, MN, MN23, MN45, u2, u23, u223, B1, B2, B23, br, br1, br2, br3, br32, rt, rt2, sign, k, r, mvc[3], a, l, alpha, beta, abc[3], m11, m21, lamb, Ptol_abs, P_frac, P_thresh, x_thresh, expect, sig, x_lo, x_hi, P, int_abserr;
    gsl_sf_result re_ln_gamma, im_ln_gamma;
    gsl_integration_workspace * int_workspace = gsl_integration_workspace_alloc(subdiv_lim);
    gsl_function p;
    
    /* Constants that depend only on M and N */
    NN = Nd*Nd;
    NN1 = Nd*(Nd+1);
    M1 = M-1;
    MN = M*Nd;
    MN23 = (MN+2)*(MN+3);
    MN45 = (MN+4)*(MN+5);
    /* Constants derived from those above, to be used later */
    u2 = 2*NN1*M*M/M1/MN23;
    br = MN*(Nd+4)-5*Nd-2;
    u23 = M1/M*MN45/4/br;
    B1 = 8/NN1/M1*MN23/MN45/MN45*br*br;
    B2 = 3/NN1/M1*MN23/MN45/(MN+6)/(MN+7)*(MN*MN*MN*(Nd+1)+MN*MN*(3*NN+68*Nd+125)-MN*(93*NN+245*Nd+32)+84*NN+48*Nd+24);
    B23 = B2+3;
    /* Even more constants derived from those above, used by Pearson I and IV */
    br2 = 2*B2-3*B1-6;
    br3 = 10*B2-12*B1-18;
    br32 = br3/br2;
    
    /* Decide which probability density function to use */
    k = B1*B23*B23/4/(4*B2-3*B1)/br2;
    
    /* Pearson IV if 0<=k<=1 */
    if (k>=0 && k<=1)
    {
      /* More constants derived from those above, dependent ultimately only on M and N */
      r = br32-2;
      mvc[0] = br32/2;
      mvc[1] = r*(2-r)*sqrt(B1/(16*(r-1)-B1*(r-2)*(r-2)));
      a = 4/sqrt(u2*(16*(r-1)-B1*(r-2)*(r-2))); /* This is 1/a in the paper*/
      l = ((r-2)/u23/4-1)*a; /* This is -l/a in the paper */
      
      /* Constant factor in probability density function, also dependent only on M and N */
      status = gsl_sf_lngamma_complex_e(mvc[0],mvc[1]/2,&re_ln_gamma,&im_ln_gamma);
      mvc[2] = 2*re_ln_gamma.val - gsl_sf_lngamma(2*mvc[0]-1) - log(2)*(2-2*mvc[0]) - log(pi); /* Don't try to take mvc[2] out of the integration as a constant factor, because you will have to use exp(mvc[2]), which will be a ridiculously small number */
      
      /* Probability density function to be integrated, giving cumulative probability */
      double pearson4 (double xla, void * mvc_ptr_void)
      {
        double *mvc_ptr = (double *)mvc_ptr_void;
        double p4 = exp(*(mvc_ptr+2) - *mvc_ptr*log(1+xla*xla) - *(mvc_ptr+1)*atan(xla));
        return p4;
      }
      /* Prepare numerical integrals */
      p.function = &pearson4;
      p.params = &mvc;
    }
    
    /* Pearson VI if k>1 */
    else if (k>1)
    {
      /* More constants derived from those above, dependent only on M and N */
      u223 = NN1*M/MN23*MN45/2/br;
      rt = 4+sqrt(16+(4+1/u2)*B1);
      alpha = u23+u223*(((u223*8-1)*u23+1)*rt+4)-1; /* alpha in MNRAS letter */
      beta = 3+2*rt/B1;
      abc[0] = alpha-1;
      abc[1] = alpha+beta;
      a = 1.0;
      l = alpha/(beta-1)-1; /* This is -delta in the paper */
      
      /* Constant factor in probability density function, also dependent only on M and N */
      abc[2] = -gsl_sf_lnbeta(alpha,beta);
      
      /* Probability density function to be integrated, giving cumulative probability */
      double pearson6 (double x, void * abc_ptr_void)
      {
        double *abc_ptr = (double *)abc_ptr_void;
        double p6 = exp(*(abc_ptr+2) + *abc_ptr*log(x) - *(abc_ptr+1)*log(x+1));
        return p6;
      }
      /* Prepare numerical integrals */
      p.function = &pearson6;
      p.params = &abc;
    }
    
    /* Pearson I if k<0*/
    else
    {
      br1 = B23/u23;
      rt2 = fabs(br3)/br3*sqrt(br1*br1+4*u2*(3*B1-4*B2)*br2);
      sign = fabs(br32)/br32;
      m11 = 1-br32/2-(1+br32/2)*sign*br1/rt2; /* This is m1+1 in the paper*/
      m21 = 2-br32-m11; /* This is m2+1 in the paper*/
      a = sign*br2/rt2;
      l = (sign*br1/rt2+1)/2;
      lamb=br1/(B1-B2+1)/3+1;
    }
    
    /* Steps common to Pearson I, IV and VI distributions */
    /* Lower threshold for P that we are aiming for (cumulative probability function) */
    P_frac = gsl_sf_erf(s/sqrt(2));
    /* Absolute tolerance on P-values*/
    Ptol_abs = Ptol*(1-P_frac)/2;
    
    /* Integrations to find upper and lower bounds for SK thresholds by trying out different P thresholds */
    for (ul=-1; ul<=1; ul+=2)
    {
      /* Threshold of cumulative probability distribution that we are aiming for */
      P_thresh = (1-ul*P_frac)/2;
      /* Expectation value of integrating variable, for Gaussian signal */
      expect = a+l;
      /* 1-sigma deviation in integrating variable, for Gaussian signal */
      sig = sqrt(u2)*fabs(a);
      /* Initial trial values of SK thresholds are symmetric thresholds */
      x_thresh = expect + ul*s*sig;
      /* Pearson I and VI functions not defined below 0, so make sure the trial threshold is not negative in these cases */
      if ((k<0 || k>1) && x_thresh<0)
        x_thresh = 0.0;
      /* Pearson I function not defined above 1, so make sure the trial threshold is not above 1 in this case */
      if (k<0 && x_thresh>1)
	x_thresh = 1.0;
      fill_lo = 0;
      fill_hi = 0;
      n = 0;
      do
      {
        if (k>=0)
	{
	  status = gsl_integration_qagiu(&p,x_thresh,Ptol_abs,0.0,subdiv_lim,int_workspace,&P,&int_abserr);
	}
	else
	{
	  P = 1-gsl_sf_beta_inc(m11,m21,x_thresh); /* This incomplete beta function is already normalised by the complete one */
	}
        if (P>P_thresh)
        {
          x_lo = x_thresh;
          x_thresh += sig; /* Steps of 1-sigma in SK threshold */
	  if (k<0 && x_thresh>1)
	    x_thresh = 1.0;
          fill_lo = 1;
        }
        else
        {
          x_hi = x_thresh;
          x_thresh -= sig;
	  if ((k<0 || k>1) && x_thresh<0)
	    x_thresh = 0.0;
          fill_hi = 1;
        }
        n++;
      }
      while ((fill_lo==0 || fill_hi==0) && n<max_its);
      /* Integrations to find SK thresholds by converging on P thresholds */
      n = 0;
      do
      {
        x_thresh = (x_lo+x_hi)/2;
        if (k>=0)
	{
	  status = gsl_integration_qagiu(&p,x_thresh,Ptol_abs,0.0,subdiv_lim,int_workspace,&P,&int_abserr);
	}
	else
	{
	  P = 1-gsl_sf_beta_inc(m11,m21,x_thresh);
	}
        if (P>P_thresh)
          x_lo = x_thresh;
        else
          x_hi = x_thresh;
        n++;
      }
      while (fabs(P-P_thresh)>Ptol_abs && n<max_its);
      sk_lims[(ul+1)/2] = (x_thresh-l)/a+lamb;
    }
    
    sk_lims_float[0] = (float)sk_lims[0];
    sk_lims_float[1] = (float)sk_lims[1];
    
    gsl_integration_workspace_free(int_workspace);
    return 0;
  }
}
