nclude <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <getopt.h>
#include <fftw3.h>


#define LIM 256

// STRUCTS
// The struct that contains the Jones matrices; frequency and the Jones matrix in complex numbers

struct jones {
  float freq;
  fftwf_complex a[4];
};



//FUNCTIONS

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

//==================================================================================
// MAIN

int main(int argc,char *argv[]) {
  
 int i,arg=0;
 int njones;
 FILE *jonesfile;
 char jonesfname[LIM],line[LIM];
 float obsfreq;
 struct jones *jm,teststruct;
 

 //getting the parameters from the command line
 while ((arg=getopt(argc,argv,"j:"))!=-1) {
    switch (arg) {
    
    case 'j':
      strcpy(jonesfname,optarg);
      break;
      
    default:
      return 0;

    }
 }
 
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
 
 //the frequency of the obs
 obsfreq=1459.;


 teststruct=jmcalculator(obsfreq,jm,njones);
 printf("tester 1 is %f\n",teststruct.freq);
 printf("tester 1 is %e\n",teststruct.a[0][0]);
 printf("tester 1 is %e\n",teststruct.a[0][1]);
 printf("tester 1 is %e\n",teststruct.a[1][0]);
 printf("tester 1 is %e\n",teststruct.a[1][1]);
 printf("tester 1 is %e\n",teststruct.a[2][0]);
 printf("tester 1 is %e\n",teststruct.a[2][1]);
 printf("tester 1 is %e\n",teststruct.a[3][0]);
 printf("tester 1 is %e\n",teststruct.a[3][1]);
 //computing the requested Jones matrix
 
 //jmatrix=jmcalculator(obsfreq,jonesarr,njones);

 //int w;
 //for(w=0;w<8;w++){printf("the element is %le\n",jmatrix[w]);} 
 // free memory blocks
// free (jonesarr);
 
 
 //close files
 
 printf("Everything works\n") ;
 
 system("PAUSE");
 return 0;   
    
}
