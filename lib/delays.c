// For these routines to work, the calc files needs to contain the reference telescope, followed by the current telescope.
// So there needs to be a calc file for every telescope.

#include "delays.h"

/* Sesame Open. Opens file with proper check. */
FILE* Sopen(char *Fin, char *how)
{
  FILE *FileIn;
  if ((FileIn  = fopen(Fin, how)) == NULL)
    {
      fprintf(stderr, "Error opening file %s\n", Fin);
      exit(0);
    }
  return(FileIn);
}

// Read Calc file
// The first time this is called, set SetCountCalcSteps to 1. This will then determine how far into the calc file to read.
// The next times this routine is called, CountCalcSteps should be increased by 1 before calling and SetCountCalcSteps should be 0
// Skipbins contains the skipbins from the input + skipbins to skip to the starttime of the last telescope + clockdelays
// Samptime is in ns
void ReadCalcfile(calc_type *Calc, double StartMJD, long SkipBins, double SampTime, int SetCountCalcSteps)
{
  int i, j, MJDint, MJDsec;
  double Refdelay, RefTelPoly[NPOLYS];
  FILE *IO;
  char line[MaxChars], Dummy[MaxChars];
  int RefTel = 0;
  int NTelescopes = 2;
  
  printf("Opening Calcfile %s\n", Calc->filename);
  IO = Sopen(Calc->filename, "r");
  for (i=0; i<18+NTelescopes; i++) fgets(line, MaxChars, IO); // Skip first lines
  fgets(line, MaxChars, IO); // Get Calc MJDint
  sscanf(line, "%18c %d\n", Dummy, &MJDint);
  fgets(line, MaxChars, IO); // Get Calc MJDsec
  sscanf(line, "%18c %d\n", Dummy, &MJDsec);
  Calc->CalcMJD = MJDint + MJDsec/86400.;
  // Calc offset (double) between start of calc and start of correlations in seconds with respect to reference telescope
  if (SetCountCalcSteps) {
    Calc->Offset = (long) ((StartMJD - Calc->CalcMJD)*86400 + 0.5) + SkipBins*SampTime/1.E9;
    Calc->CountCalcSteps = (int)(Calc->Offset/CALCSTEPSIZE);
  }
  // Read polynomials. First for reftel.
  for (i=0; i<(2*6*NTelescopes+2)*Calc->CountCalcSteps; i++) fgets(line, MaxChars, IO); // Skip to right timesection
  fgets(line, MaxChars, IO);
  sscanf(line, "%23c %lf %lf %lf %lf %lf %lf", Dummy, &RefTelPoly[0], &RefTelPoly[1], &RefTelPoly[2], &RefTelPoly[3], &RefTelPoly[4], &RefTelPoly[5]);
  //  printf("%e %e %e %e %e %e\n", RefTelPoly[0], RefTelPoly[1], RefTelPoly[2], RefTelPoly[3], RefTelPoly[4], RefTelPoly[5]);
  for (j=0; j<5; j++) fgets(line, MaxChars, IO);
  fgets(line, MaxChars, IO);
  sscanf(line, "%23c %lf %lf %lf %lf %lf %lf", Dummy, &Calc->Poly[0], &Calc->Poly[1], &Calc->Poly[2], &Calc->Poly[3],
	 &Calc->Poly[4], &Calc->Poly[5]);
  //  printf("%e %e %e %e %e %e\n", Calc->Poly[0], Calc->Poly[1], Calc->Poly[2], Calc->Poly[3], Calc->Poly[4], Calc->Poly[5]);
  for (j=0; j<5; j++) fgets(line, MaxChars, IO);
  fclose(IO);
  
  for (j=0; j<NPOLYS; j++)
    Calc->Poly[j] -= RefTelPoly[j];
}

// Read the GPS files
void ReadGPSFiles(gps_type *GPS, double StartMJD)
{
  int i, Done;
  FILE *IO;
  char line[MaxChars];
  float MJD, ClockOffset_read;
  float MJD_prev=0, ClockOffset_read_prev;
  float ClockOffset_Ref;
  double ClockDrift_Ref;
  Done = 0;
  GPS->ClockDrift=0;
  printf("Opening GPS file %s\n", GPS->filename);
  IO = Sopen(GPS->filename, "r");
  while (fgets(line, MaxChars, IO) && !Done) {
    if (line[0] != '#') {
      sscanf(line, "%f %f", &MJD, &ClockOffset_read);
      if (StartMJD < MJD) {
	GPS->ClockOffset = (ClockOffset_read - ClockOffset_read_prev) * (StartMJD-MJD_prev)/(MJD-MJD_prev) + ClockOffset_read_prev;
	GPS->ClockDrift += (ClockOffset_read - ClockOffset_read_prev)/86400; // in us / s
	Done = 1;
      }
      MJD_prev = MJD;
      ClockOffset_read_prev = ClockOffset_read;
    }
  }
  if (!Done) {fprintf(stderr, "Can't get clockoffset!\n"); GPS->ClockOffset = 0;}
  if (GPS->ClockOffset == 0) {
    fprintf(stderr, "Could not find Clockoffset with MJD=%f in %s\n", StartMJD, GPS->filename);
    exit(0);
  }
  // Read the reference clock.
  Done = 0;
  ClockDrift_Ref=0;
  IO = Sopen(GPS->filenameref, "r");
  while (fgets(line, MaxChars, IO) && !Done) {
    if (line[0] != '#') {
      sscanf(line, "%f %f", &MJD, &ClockOffset_read);
      if (StartMJD < MJD) {
	ClockOffset_Ref = (ClockOffset_read-ClockOffset_read_prev) * (StartMJD-MJD_prev)/(MJD-MJD_prev) + ClockOffset_read_prev;
	ClockDrift_Ref += (ClockOffset_read - ClockOffset_read_prev)/86400; // in us / s
	Done = 1;
      }
      MJD_prev = MJD;
      ClockOffset_read_prev = ClockOffset_read;
    }
  }
  if (!Done) {fprintf(stderr, "Can't get clockoffset!\n"); ClockOffset_Ref = 0;}
  if (ClockOffset_Ref == 0) {
    fprintf(stderr, "Could not find Clockoffset with MJD=%f in %s\n", StartMJD, GPS->filenameref);
    exit(0);
  }
  // Move reference of clockoffset and clockdrift to telescope 0, so that telescope 0 has no offset and no drift
  GPS->ClockOffset -= ClockOffset_Ref;
  GPS->ClockDrift -= ClockDrift_Ref;
}
