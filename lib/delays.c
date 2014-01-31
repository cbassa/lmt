// For these routines to work, the calc files needs to contain the reference telescope, followed by the current telescope.
// So there needs to be a calc file for every telescope.

// Read Calc file
// The first time this is called, set SetCountCalcSteps to 1. This will then determine how far into the calc file to read.
// The next times this routine is called, CountCalcSteps should be increased by 1 before calling and SetCountCalcSteps should be 0
void ReadCalcfile(Calc_type *Calc, double StartMJD, long SkipBins, double SampTime, int SetCountCalcSteps)
{
  int i, j, MJDint, MJDsec;
  double Refdelay, RefTelPoly[6];
  FILE *IO;
  char line[MaxChars], Dummy[MaxChars];
  int RefTel = Fileinfo->RefTel;
  int NTelescopes = 2;
  
  IO = Sopen(Calc->filename, "r");
  for (i=0; i<18+Fileinfo->NTelescopes; i++) fgets(line, MaxChars, IO); // Skip first lines
  fgets(line, MaxChars, IO); // Get Calc MJDint
  sscanf(line, "%18c %d\n", Dummy, &MJDint);
  fgets(line, MaxChars, IO); // Get Calc MJDsec
  sscanf(line, "%18c %d\n", Dummy, &MJDsec);
  Calc->CalcMJD = MJDint + MJDsec/86400.;
  // Calc offset (double) between start of calc and start of correlations in seconds with respect to reference telescope
  if (SetCountCalcSteps) {
    Calc->Offset = (long) ((Fileinfo->StartMJD[RefTel] - Calc->CalcMJD)*86400 + 0.5) + Fileinfo->SkipBins[RefTel]*Fileinfo->SampTime[RefTel]/1.E9;
    Calc->CountCalcSteps = (int)(Calc->Offset/CALCSTEPSIZE);
  }
  // Read polynomials. First for reftel.
  for (i=0; i<(2*6*Fileinfo->NTelescopes+2)*Calc->CountCalcSteps; i++) fgets(line, MaxChars, IO); // Skip to right timesection
  fgets(line, MaxChars, IO);
  sscanf(line, "%23c %lf %lf %lf %lf %lf %lf", Dummy, &RefTelPoly[0], &RefTelPoly[1], &RefTelPoly[2], &RefTelPoly[3], &RefTelPoly[4], &RefTelPoly[5]);
  for (j=0; j<5; j++) fgets(line, MaxChars, IO);
  fgets(line, MaxChars, IO);
  sscanf(line, "%23c %lf %lf %lf %lf %lf %lf", Dummy, &Calc->Poly[i][0], &Calc->Poly[i][1], &Calc->Poly[i][2], &Calc->Poly[i][3],
	 &Calc->Poly[i][4], &Calc->Poly[i][5]);
  for (j=0; j<5; j++) fgets(line, MaxChars, IO);
  fclose(IO);
  
  for (j=0; j<NPolys; j++)
    Calc->Poly[j] -= RefTelPoly[j];
  free(RefTelPoly);
}

// Read the GPS files
void ReadGPSFiles(Fileinfo_type *Fileinfo)
{
  int i, Done;
  FILE *IO;
  char line[MaxChars];
  float MJD, ClockOffset_read;
  float MJD_prev=0, ClockOffset_read_prev;
  for (i=0; i<Fileinfo->NTelescopes; i++) {
    Done = 0;
    IO = Sopen(Fileinfo->GPSFiles[i], "r");
    while (fgets(line, MaxChars, IO) && !Done) {
      if (line[0] != '#') {
        sscanf(line, "%f %f", &MJD, &ClockOffset_read);
        if (Fileinfo->StartMJD[i] < MJD) {
          Fileinfo->ClockOffset[i] = (ClockOffset_read-ClockOffset_read_prev) * (Fileinfo->StartMJD[i]-MJD_prev)/(MJD-MJD_prev) + ClockOffset_read_prev;
          Fileinfo->ClockDrift[i] += (ClockOffset_read - ClockOffset_read_prev)/86400; // in us / s
	  Done = 1;
        }
        MJD_prev = MJD;
        ClockOffset_read_prev = ClockOffset_read;
      }
    }
    if (!Done) {fprintf(stderr, "Can't get clockoffset!\n"); Fileinfo->ClockOffset[i] = 0;}
    if (Fileinfo->ClockOffset[i] == 0) {
      fprintf(stderr, "Could not find Clockoffset with MJD=%f in %s\n", Fileinfo->StartMJD, Fileinfo->GPSFiles[0]);
      exit(0);
    }
  }
  // Move reference of clockoffset and clockdrift to telescope 0, so that telescope 0 has no offset and no drift
  for (i=1; i<Fileinfo->NTelescopes; i++) {
    Fileinfo->ClockOffset[i] -= Fileinfo->ClockOffset[0];
    Fileinfo->ClockDrift[i] -= Fileinfo->ClockDrift[0];
  }
  Fileinfo->ClockOffset[0] = 0;
  Fileinfo->ClockDrift[0] = 0;
}

