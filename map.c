#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <getopt.h>
#include "filterbank.h"


#define LIM 256


// STRUCTS
// The struct that contains the Jones matrices; frequency and the Jones matrix in complex numbers

struct jones {
  float freq;
  fftwf_complex a[4];
};



//FUNCTIONS

// Read a line of maximum length int lim from file FILE into string s
int fgetline(FILE *file,char *s,int lim)
{
  int c,i=0;

  while (--lim > 0 && (c=fgetc(file)) != EOF && c != '\n')
    s[i++] = c;
  if (c == '\t')
    c=' ';
  if (c == '\n')
    s[i++] = c;
  s[i] = '\0';
  return i;
}

//The function that computes the Jones matrix for a specific frequency. Its input are 
//"frequency, the Jones matrix array from the input file, the number of frequencies it includes"
struct jones jmcalculator (float freq, struct jones *fjones, int njones){
       
 int counter=0,i,j; // counter is the counter that will give the frequency just above obsfreq
 struct jones jonesm;
       
 if(freq<fjones[0].freq)
    {
    fprintf(stderr,"The observation frequency is lower than the minimum frequency in the Jones matrix... Exiting\n");
    exit;                      
    }  
  else if(freq>fjones[njones-1].freq)
    {
    fprintf(stderr,"The observation frequency is higher than the maximum frequency in the Jones matrix... Exiting\n");
    exit;                      
    }  
  else
    {
    while(fjones[counter].freq<freq){counter++;}                           
    }

 //computing the Jones matrix, jonesm, that will be used for obsfreq. Usage of linear interpolation
 jonesm.freq=freq;
 
 for(i=0;i<4;i++) {
        for(j=0;j<2;j++)
        jonesm.a[i][j]=((fjones[counter].a[i][j]-fjones[counter-1].a[i][j])/(fjones[counter].freq-fjones[counter-1].freq))*(freq-fjones[counter-1].freq)+fjones[counter-1].a[i][j];                  
 }
 
    return jonesm;
}


// MAIN


int main(int argc,char *argv[])
{
  int arg=0;
  int i,j,bytes_read,njones;
  fftwf_complex *rp1,*rp2,*cp1,*cp2;
  FILE *infile,*outfile,*jonesfile;
  char infname[LIM],outfname[LIM],jonesfname[LIM],line[LIM];
  struct filterbank fbin,fbout;
  struct jones *jm,*fjm;
  float obsfreq;
  

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
  printf("map: applying jones matrix %s\n",jonesfname);

  // Write header struct
  fwrite(&fbout,1,sizeof(struct filterbank),outfile);

  // Allocate buffers
  rp1=fftwf_malloc(sizeof(fftwf_complex)*fbin.nchan);
  rp2=fftwf_malloc(sizeof(fftwf_complex)*fbin.nchan);
  cp1=fftwf_malloc(sizeof(fftwf_complex)*fbout.nchan);
  cp2=fftwf_malloc(sizeof(fftwf_complex)*fbout.nchan);
  
  
  //This section of the code reads the Jones matrix file and computes the 
  //Jones matrix specifically for the frequency of the obs file
  
  //Open the file with the Jones matrix
  jonesfile=fopen(jonesfname,"r");
  
  // Check if input file exists
  if (jonesfile==NULL) {
    fprintf(stderr,"Error opening %s\n",jonesfname);
    exit;
 }  
  
  // Count elements in polcal file
  i=0;
  while (fgetline(jonesfile,line,LIM)>0) 
    i++;
  rewind(jonesfile);
  njones=i;  
   
  // Allocate the memory for all the jones matrices retrieved from the file
  jm=(struct jones *) malloc(sizeof(struct jones)*njones); 
  
  // Read jones matrices
  for (i=0;i<njones;i++) {
    fgetline(jonesfile,line,LIM);
    sscanf(line,"%f %e %e %e %e %e %e %e %e",&jm[i].freq,
	   &jm[i].a[0][0],
	   &jm[i].a[0][1],
	   &jm[i].a[1][0],
	   &jm[i].a[1][1],
	   &jm[i].a[2][0],
	   &jm[i].a[2][1],
	   &jm[i].a[3][0],
	   &jm[i].a[3][1]);
  }  
  
  //Closing the jonesfile
  fclose(jonesfile);
  
  // Creating an array of jones structs that will contain the Jones matrices for all frequency bins
  fjm=(struct jones *) malloc(sizeof(struct jones)*fbin.nchan); 
  
  for(i=0;i<fbin.nchan;i++) {
      obsfreq=fbin.freq+fbin.bw*(float)i/(float)(fbin.nchan-1)-0.5*fbin.bw;
      fjm[i]=jmcalculator(obsfreq,jm,njones);                          
  }
  printf("Applying Jones matrices in frequency channels...\n");
        
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
      cp1[j][0]=fjm[j].a[0][0]*rp1[j][0]-fjm[j].a[0][1]*rp1[j][1]+fjm[j].a[1][0]*rp2[j][0]-fjm[j].a[1][1]*rp2[j][1];
      cp1[j][1]=fjm[j].a[0][0]*rp1[j][1]+fjm[j].a[0][1]*rp1[j][0]+fjm[j].a[1][0]*rp2[j][1]+fjm[j].a[1][1]*rp2[j][0];
      cp2[j][0]=fjm[j].a[2][0]*rp1[j][0]-fjm[j].a[2][1]*rp1[j][1]+fjm[j].a[3][0]*rp2[j][0]-fjm[j].a[3][1]*rp2[j][1];
      cp2[j][1]=fjm[j].a[2][0]*rp1[j][1]+fjm[j].a[2][1]*rp1[j][0]+fjm[j].a[3][0]*rp2[j][1]+fjm[j].a[3][1]*rp2[j][0];
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
  free(jm);
  free(fjm);

  return 0;
}
