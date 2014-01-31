#ifndef CALC_H
#define CALC_H

#define MaxChars 999

// Define Calc structure
typedef struct {
  double **Poly;
  char *filename;
  double Offset; // Offset from Calc to start of correlation in seconds
  int CountCalcSteps; // Count the number of CALCSTEPSIZE into the observation
  double CalcMJD; // The first MJD in the Calc file
} calc_type;

void ReadCalcfile(Calc_type *Calc, Fileinfo_type *Fileinfo, int SetCountCalcSteps);

#endif
