#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <getopt.h>
#include "filterbank.h"


#define LIM 256


// STRUCTS


//FUNCTIONS


// MAIN


int main(int argc,char *argv[]) {

 int arg=0,bytes_read;
 FILE *infileT1,*infileT2,*outfile;
 char infT1name[LIM],infT2name[LIM],outfname;
 struct filterbank fbinT1,fbinT2,fbout;
 fftwf_complex *rp1T1,*rp2T1,*rp1T2,*rp2T2,*cp1,*cp2,*cp3,*cp4,*cp5,*cp6,*cp7,*cp8,;
 
 
// Getting the parameters from the command line
  while ((arg=getopt(argc,argv,"i:j:o:"))!=-1) {
    switch (arg) {

    case 'i':
      strcpy(infT1name,optarg);
      break;

    case 'j':
      strcpy(infT2name,optarg);
      break;
      
    case 'o':
      strcpy(outfname,optarg);
      break;

    default:
      return 0;

    }
  }
    
   // Open input files
  infileT1=fopen(infT1name,"w");
  infileT2=fopen(infT2name,"w");
  
  // Check if input files exist
  if (infileT1==NULL) {
    fprintf(stderr,"Error opening %s\n",infT1name);
    exit;
  }
  
   if (infileT2==NULL) {
    fprintf(stderr,"Error opening %s\n",infT2name);
    exit;
  }

  // Open output file
  outfile=fopen(outfname,"w");
  
  // Check if output file exists
  if (outfile==NULL) {
    fprintf(stderr,"Error opening %s\n",outfname);
    exit;
  }
   
   // Read headers from the two files
  bytes_readT1=fread(&fbinT1,1,sizeof(struct filterbank),infT1name);
  bytes_readT2=fread(&fbinT2,1,sizeof(struct filterbank),infT2name);

  // Copy filterbank struct (we must create from scratch I think, at least partially.not placing anything here for the moment)
  
  
 
  //Allocating buffers
  rp1T1=fftwf_malloc(sizeof(fftwf_complex)*fbinT1.nchan);
  rp2T1=fftwf_malloc(sizeof(fftwf_complex)*fbinT1.nchan);
  rp1T2=fftwf_malloc(sizeof(fftwf_complex)*fbinT2.nchan);
  rp2T2=fftwf_malloc(sizeof(fftwf_complex)*fbinT2.nchan);
  
  
  // Print information
  printf("Correlator started...\n");
  
  
  //Reading the data from the input files
  for (;;) {
    // Read buffers
    bytes_readT1=fread(rp1T1,sizeof(fftwf_complex),fbinT1.nchan,infileT1);
    bytes_readT1=fread(rp2T1,sizeof(fftwf_complex),fbinT1.nchan,infileT1);
    bytes_readT2=fread(rp1T2,sizeof(fftwf_complex),fbinT2.nchan,infileT2);
    bytes_readT2=fread(rp2T2,sizeof(fftwf_complex),fbinT2.nchan,infileT2);
   
   
   
   
   return 0;
}
