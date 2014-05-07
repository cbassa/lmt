/* This is a new version of ppredict.c, made by Ben Stappers
   This new version is adapted for use with add_phased

   Roy Smits 17 Nov 2011
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h>

#include "ppolyco.c"
#include "mpolyco.c"
#ifdef __linux__
#include <getopt.h>
#endif

#define PROG "PPREDICT:    "

void usage(int status)
{
  printf("%s A program to predict the period for a particular pulsar at a certain time and site\n");
  printf("------------------------------------------------------------------------------------\n");
  printf("%sUsage:    ppredict ([-i intmjd -f fracmjd] || [-u UT time]) [options] -s site psrname\n",PROG);
  printf("------------------------------------------------------------------------------------\n");
  printf("%sOptions:  -h prints this help information\n",PROG);
  printf("%s          -c Give the number of coefficients to use for the polyco, i.e. 12 (default)\n",PROG);    
  printf("%s          -f Give the fractional part of the MJD\n",PROG);
  printf("%s          -i Give the integer part of the MJD\n",PROG);
  printf("%s          -n Give the span for each polyco, i.e. 180s (this is the default, should be an integer)\n",PROG);
  printf("%s          -s Give the site for which the period is to be calculated (default is WSRT)\n",PROG);
  printf("%s          -u Give the UT time (yyyy:mm:dd:hh:mm:ss.sss). Mutually exclusive with -i & -f options\n",PROG);
  printf("%s          -v Verbose mode\n",PROG);
}

int verb1,verb2,verb3,verb4;

void predict(int intmjd, double fracmjd, char *psr, char *parfname, char site[2], double *phase, double *period, double *dm)
{
  int    coeff_set,day,err,hour,i,min,month,ncoeff,nlen,nspan,site_set,span_set,ut_set,year;
  float  sec;
  double frac,midpobs,mjd,refph;
  char   c,*unfname,utchar[32], *psrname;
  
  /* Variable Initialization */
  ut_set = site_set = coeff_set = span_set = 0;
  nspan = 0;
  ncoeff = 0;
  midpobs = 0.0;
  refph = 0.0;
  *dm = 0.0;

  verb1 = verb2 = verb3 = 0;

  /* Converting the pulsar name to filenames*/
  nlen = strlen(psr);
  unfname = (char *) calloc(nlen,1);
  psrname = (char *) calloc(nlen,1);
  strcpy(unfname,psr);
  strcpy(psrname,psr);
  //  printf("# PSRname: %s\n", psrname);
  mjd = (double) (intmjd);
  sla_djcl_(&mjd,&year,&month,&day,&frac,&err);
  fracmjd += frac;
  if (fracmjd >= 1.0) {fracmjd -= 1.0; day += 1;}
  hour = floor(24 * fracmjd);
  min  = floor(60 * ((fracmjd * 24) - hour));
  sec  = floor(60 * ((60 * ((24 * fracmjd) - hour)) - min));
  // Use tempo2 to make a polyco based on the par file given by parfname
  mpolyco_t1(unfname,psrname,intmjd,fracmjd,site,&nspan,&ncoeff,parfname);
  // Calculate pulsar period and phase from the polyco
  ppolyco(unfname,intmjd,fracmjd,&midpobs,&refph,dm);
  mjd = intmjd + fracmjd;
  *phase = refph-(long)refph;
  *period = midpobs;
}
