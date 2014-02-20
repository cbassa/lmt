#ifndef CALC_H
#define CALC_H

#include <stdio.h>
#include <stdlib.h>

#define MaxChars 999
#define CALCSTEPSIZE 120
#define NPOLYS 6 // When changing this value, check the code (or better yet, just don't change this value)

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

void ReadCalcfile(calc_type *Calc, double StartMJD, long SkipBins, double SampTime, int SetCountCalcSteps);

#endif
