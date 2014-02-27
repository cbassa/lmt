#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <getopt.h>


#define LIM 256

//function prototypes
double *jmcalculator (double freq, double *array, int freqsno);

int main(int argc,char *argv[])
{  
 int arg=0;
 int count1,jmfnofreq;
 FILE *jonesfile;
 char jonesfname[LIM];
 double *jonesarr,temp,obsfreq,*jmatrix;
 

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
 
 //Reading the jones matrix file
 count1=0;   //count1 is the counter of the Jones matrix file elements
 jonesarr=malloc(sizeof(double));
 while(fscanf(jonesfile,"%le",&temp)!=EOF)
    {         
    jonesarr=realloc(jonesarr,(count1+1)*sizeof(double));
    jonesarr[count1]=temp;
    //printf("it is %.7le\n",jonesarr[count1]);
    count1++;
    }
 
 //the number of frequencies in the Jones matrix file
 jmfnofreq=count1/9;
  
 //the frequency of the obs
 obsfreq=1459.;

 //computing the requested Jones matrix
 jmatrix=malloc(7*sizeof(double));
 jmatrix=jmcalculator(obsfreq,jonesarr,jmfnofreq);
 
 // free memory blocks
 free (jonesarr);
 free (jmatrix);
 
 //close files
 close (jonesfile);
 
 system("PAUSE");
 return 0;   
    
}

//FUNCTIONS

//the function that computes the Jones matrix for a specific frequency. Its input are 
//"frequency, the Jones matrix array from the input file, the number of frequencies it includes"
double *jmcalculator (double freq, double *array, int freqsno){
       
 int countfun1=0,i; // count2 is the counter that will give the frequency just above obsfreq
 double freqb1,freqb2;
 static double jonesm[7];
       
 if(freq<array[0])
    {
    fprintf(stderr,"The observation frequency is lower than the minimum frequency in the Jones matrix... Exiting\n");
    exit;                      
    }  
  else if(freq>array[(freqsno-1)*9])
    {
    fprintf(stderr,"The observation frequency is higher than the maximum frequency in the Jones matrix... Exiting\n");
    exit;                      
    }  
  else
    {
    while(array[countfun1*9]<freq){countfun1++;}                           
    }
    
 freqb1=array[9*(countfun1-1)];
 freqb2=array[9*countfun1];

 //computing the Jones matrix, jonesm, that will be used for obsfreq. Usage of linear interpolation
 for(i=0;i<8;i++)
    {
    jonesm[i]=((array[countfun1*9+1+i]-array[(countfun1-1)*9+1+i])/(freqb2-freqb1))*(freq-freqb1)+array[(countfun1-1)*9+1+i]; 
    }
    return jonesm;
}
