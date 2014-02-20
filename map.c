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
  int i,j,bytes_read,count1,count2,countf,jmfnofreq;
  fftwf_complex *rp1,*rp2,*cp1,*cp2;
  FILE *infile,*outfile,*jonesfile;
  char infname[LIM],outfname[LIM],jonesfname[LIM];
  double *jonesarr,temp,obsfreq,freqb1,freqb2,*jmffreqarr,jonesm[7];
  struct filterbank fbin,fbout;

  // Getting the parameters from the command line
  while ((arg=getopt(argc,argv,"i:o:j:"))!=-1) {
    switch (arg) {

    case 'i':
      strcpy(infname,optarg);
      break;

    case 'o':
      strcpy(outfname,optarg);
      break;
      
    case 'j':
      strcpy(jonesfname,optarg);
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
  
  
  //============================================================================
  //This section of the code reads the Jones matrix file and computes the 
  //Jones matrix specifically for the frequency of the obs file
  
  //Open the file with the Jones matrix
 jonesfile=fopen(jonesfname,"r");
 
 // Check if input file exists
 if (jonesfile==NULL) {
    fprintf(stderr,"Error opening %s\n",jonesfname);
    exit;
 }
 
 //Reading the jones matrix file
 count1=0;   //count1 is the counter of the Jones matrix file elements
 jonesarr=malloc(sizeof(double));
 while(fscanf(jonesfile,"%le",&temp)!=EOF)
    {         
    jonesarr=realloc(jonesarr,(count1+1)*sizeof(double));
    jonesarr[count1]=temp;
    count1++;
    }
 
 //the number of frequencies in the Jones matrix file
 jmfnofreq=count1/9;
  
 //creating an array with the frequencies in the jones matrix file
 countf=0;
 jmffreqarr=malloc(jmfnofreq*sizeof(double));
 while(countf<jmfnofreq)
    {
    jmffreqarr[countf]=jonesarr[countf*9];               
    countf++;
    }
 
 //the frequency of the obs
 obsfreq=fbin.freq;

 //finding the frequencies [freqb1,freqb2] from the jmffreqarr that surround obsfreq
 count2=0;  // count2 is the counter that will give the frequencyjust above obsfreq
 if(obsfreq<jmffreqarr[0])
    {
    fprintf(stderr,"The observation frequency is lower than the minimum frequency in the Jones matrix... Exiting\n");
    exit;                      
    }  
 else if(obsfreq>jmffreqarr[countf-1])
    {
    fprintf(stderr,"The observation frequency is higher than the maximum frequency in the Jones matrix... Exiting\n");
    exit;                      
    }  
 else
    {
    while(jmffreqarr[count2]<obsfreq){count2++;}                           
    }
    
 freqb1=jmffreqarr[count2-1];
 freqb2=jmffreqarr[count2];
 
 //computing the Jones matrix, jonesm, that will be used for obsfreq. Usage of linear interpolation
 i=0;
 for(i=0;i<8;i++)
    {
    jonesm[i]=((jonesarr[count2*9+1+i]-jonesarr[(count2-1)*9+1+i])/(freqb2-freqb1))*(obsfreq-freqb1)+jonesarr[(count2-1)*9+1+i];     
    }
  
  
  //============================================================================
  
  
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
      cp1[j][0]=jonesm[0]*rp1[j][0]-jonesm[1]*rp1[j][1]+jonesm[2]*rp2[j][0]-jonesm[3]*rp2[j][1];
      cp1[j][1]=jonesm[1]*rp1[j][0]+jonesm[0]*rp1[j][1]+jonesm[2]*rp2[j][1]+jonesm[3]*rp2[j][0];
      cp2[j][0]=jonesm[4]*rp1[j][0]-jonesm[5]*rp1[j][1]+jonesm[6]*rp2[j][0]-jonesm[7]*rp2[j][1];
      cp2[j][1]=jonesm[4]*rp1[j][1]+jonesm[5]*rp1[j][0]+jonesm[6]*rp2[j][1]+jonesm[7]*rp2[j][0];
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
