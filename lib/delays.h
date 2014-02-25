#ifndef CALC_H
#define CALC_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MaxChars 999
#define CALCSTEPSIZE 120
#define NPOLYS 6 // When changing this value, check the code (or better yet, just don't change this value)
#define MAXDELAYOFFSET 0.03 // This is the number of seconds to add to the beginning of each stream so that we can correct for geometric delay
#define DELTA 0.0000001

const float PI = M_PI;


// Define Calc structure
typedef struct {
  double *Poly;
  char *filename;
  double Offset; // Offset from Calc to start of correlation in seconds
  int CountCalcSteps; // Count the number of CALCSTEPSIZE into the observation
  double CalcMJD; // The first MJD in the Calc file
} calc_type;

typedef struct {
  char *filename;
  char *filenameref;
  float ClockOffset;
  double ClockDrift;
} gps_type;

typedef struct{
  long skipbins;
  long binsshifted;
  long geoshifted;
  float frac_delay_pol1;
  float frac_delay_pol2;
  float phase_offset_pol1;
  float phase_offset_pol2;
} delays_type;

FILE* Sopen(char *Fin, char *how);
void ReadCalcfile(calc_type *Calc, double StartMJD, long SkipBins, double SampTime, int SetCountCalcSteps);
void ReadGPSFiles(gps_type *GPS, double StartMJD);
void RotateVec(float *re, float *im, float angle);
double GetStartGeoOffset(calc_type *Calc, double samptime, long reducedbinnumber);
double GetTimedelay(delays_type *delays, calc_type *Calc, double samptime, long reducedbinnumber);
//double GetFringeRotation(int N, int Telescope, int freqcounter, double timedelay, Calc_type *Calc, Fileinfo_type *Fileinfo, double time, int Pol);
void Add_Geo_Delays(delays_type *delays, double startgeooffset, double samptime);
void Add_Clock_Delays(delays_type *delays, gps_type *gps, double samptime);

#endif
