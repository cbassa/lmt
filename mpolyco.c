#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>

#define MPROG "MPOLYCO: "

int verb1,verb2,verb3,verb4;


void mpolyco(fname,psrname,imjd,fmjd,nsite,nspan,ncoeff)
char  *fname,*psrname,*nsite;
int   *nspan,*ncoeff,imjd;
double fmjd;
{
  FILE     *tzfile,*datefile;

  char     unique[32],tzin[16],date[16],mvpoly[2000],genpoly[2000];

  int      maxha; 
  float    freq;
  
  /* Get these from catalogue at a later date */
  *ncoeff = 15; 
  *nspan = 120;
  maxha = 12;
  freq = 1408.0;
 
  if (verb2) printf("mp:psrname = %s\n",psrname);
  /* First we need to make the tz.in file */
  
  strcpy(tzin,psrname);
  strcat(tzin,".tz");
  tzfile = fopen(tzin,"w");

  fprintf(tzfile,"    %s   12 180  15 1408      (Defaults for nsite, maxha, nspan, ncoeff, freq)\n",nsite);
  fprintf(tzfile,"  Name     Nspan  Ncoeffs Maxha Freq (For each PSR you can override the defaults)\n");
  fprintf(tzfile,"--------------------------------------------------------------------------------\n");
  fprintf(tzfile,"%s %d %d %d %f\n",psrname,*nspan,*ncoeff,maxha,freq); 
  fclose(tzfile);

  /* Okay now we need a file to put the date in */
  
  if (verb2) printf("fname = %s\n",fname);

  strcpy(date,psrname);
  strcat(date,".dt");
  datefile = fopen(date,"w");
  fprintf(datefile,"%lf %lf\n",imjd+fmjd-0.5,imjd+fmjd+0.5);
  fclose(datefile);

  /* Actually generate the polyco */
  sprintf(genpoly,"tempo -z %s < %s >> %s",tzin,date,fname);
  //  sprintf(genpoly,"tempo -z %s < %s >> /tmp/%s",tzin,date,fname);
  if (verb2) printf("%s %s\n",MPROG,genpoly);
  printf("Generating polyco with Tempo.\n");
  system(genpoly);

  /* Move the polyco to a unique name */

  sprintf(mvpoly,"mv polyco.dat %s.polyco\n",fname);
  if (verb2) printf("%s %s\n",MPROG,mvpoly);
  system(mvpoly);
}

void mpolyco_t1(fname,psrname,imjd,fmjd,nsite,nspan,ncoeff,parfname)
char  *fname,*psrname,*nsite,*parfname;
int   *nspan,*ncoeff,imjd;
double fmjd;
{
  FILE     *tzfile,*datefile;

  char     unique[32],tzin[16],date[16],mvpoly[2000],genpoly[2000];

  int      maxha; 
  float    freq;
  
  /* Get these from catalogue at a later date */
  *ncoeff = 15; 
  *nspan = 120;
  maxha = 12;
  freq = 1408.0;
 
  if (verb2) printf("mp:psrname = %s\n",psrname);
  /* First we need to make the tz.in file */
  
  strcpy(tzin,psrname);
  strcat(tzin,".tz");
  tzfile = fopen(tzin,"w");

  fprintf(tzfile,"    %s   12 180  15 1408      (Defaults for nsite, maxha, nspan, ncoeff, freq)\n",nsite);
  fprintf(tzfile,"  Name     Nspan  Ncoeffs Maxha Freq (For each PSR you can override the defaults)\n");
  fprintf(tzfile,"--------------------------------------------------------------------------------\n");
  fprintf(tzfile,"%s %d %d %d %f\n",psrname,*nspan,*ncoeff,maxha,freq); 
  fclose(tzfile);

  /* Okay now we need a file to put the date in */
  
  if (verb2) printf("fname = %s\n",fname);

  strcpy(date,psrname);
  strcat(date,".dt");
  datefile = fopen(date,"w");
  fprintf(datefile,"%lf %lf\n",imjd+fmjd-0.5,imjd+fmjd+0.5);
  fclose(datefile);

  /* Actually generate the polyco */
  sprintf(genpoly,"tempo -f %s -z %s < %s >> %s",parfname, tzin,date,fname);
  //  sprintf(genpoly,"tempo -f %s -z %s < %s >> /tmp/%s",parfname, tzin,date,fname);
  if (verb2) printf("%s %s\n",MPROG,genpoly);
  printf("Generating polyco with Tempo.\n");
  system(genpoly);

  /* Move the polyco to a unique name */

  sprintf(mvpoly,"mv polyco.dat %s.polyco\n",fname);
  if (verb2) printf("%s %s\n",MPROG,mvpoly);
  system(mvpoly);
}

void mpolyco_t2(fname,psrname,imjd,fmjd,nsite,nspan,ncoeff,parfname)
char  *fname,*psrname,*nsite,*parfname;
int   *nspan,*ncoeff,imjd;
double fmjd;
{
  FILE     *tzfile,*datefile;

  char     mvpoly[2000],genpoly[2000];

  int      maxha; 
  float    freq;
  
  /* Get these from catalogue at a later date */
  *ncoeff = 15; 
  *nspan = 120;
  maxha = 12;
  freq = 1408.0;
 
  if (verb2) printf("mp:psrname = %s\n",psrname);

  /* Actually generate the polyco */
  sprintf(genpoly,"tempo2 -f %s -polyco \"%f %f %d %d %d %s %f\" -tempo1", parfname,imjd+fmjd-0.5,imjd+fmjd+0.5,*nspan,*ncoeff,maxha,nsite, freq);
  if (verb2) printf("%s %s\n",MPROG,genpoly);
  printf("Generating polyco with Tempo2.\n");
  system(genpoly);

  /* Move the polyco to a unique name */

  sprintf(mvpoly,"mv polyco_new.dat %s.polyco",fname);
  /*  if(!verb2) strcat(mvpoly, " > /dev/null"); */
  if (verb2) printf("%s %s\n",MPROG,mvpoly);
  system(mvpoly);
}






void mpolyco_t0(fname,psrname,imjd,fmjd,nsite,nspan,ncoeff,parfname)
char  *fname,*psrname,*nsite,*parfname;
int   *nspan,*ncoeff,imjd;
double fmjd;
{
  FILE     *tzfile,*datefile;

  char     mvpoly[2000],genpoly[2000];

  int      maxha; 
  float    freq;
  
  /* Get these from catalogue at a later date */
  *ncoeff = 15; 
  *nspan = 120;
  maxha = 12;
  freq = 1408.0;
 
  if (verb2) printf("mp:psrname = %s\n",psrname);

  /* Actually generate the polyco */
  
  sprintf(genpoly,"polyco -f %s %f %f %d %d %d %s %f", parfname,imjd+fmjd-0.5,imjd+fmjd+0.5,*nspan,*ncoeff,maxha,nsite, freq);
  if (verb2) printf("%s %s\n",MPROG,genpoly);
  printf("Generating polyco with programme 'Polyco'.\n");
  system(genpoly);

  /* Move the polyco to a unique name */

  sprintf(mvpoly,"mv polyco.dat %s.polyco",fname);
  /*  if(!verb2) strcat(mvpoly, " > /dev/null"); */
  if (verb2) printf("%s %s\n",MPROG,mvpoly);
  system(mvpoly);
}






