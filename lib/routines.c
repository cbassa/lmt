void Uppercase(char *string)
{
  int i;
  for (i=0; string[i]; i++)
    string[i] = toupper(string[i]);
}


// Convert Telescope 2-letter name to one char for tempo
// Telescope name is already converted to uppercase.
char* Telescope2Char(char *Telescope)
{
  if (!strcmp(Telescope, "EB")) return "g";
  if (!strcmp(Telescope, "WB")) return "i";
  if (!strcmp(Telescope, "JB")) return "8";
  if (!strcmp(Telescope, "NC")) return "f";
  fprintf(stderr, "Error, Telescope %s not known\n", Telescope);
  exit(0);
}

// angle in radians
void RotateVec(float *re, float *im, float angle)
{
  float x, y, cosangle, sinangle;
  x = *re;
  y = *im;
  //  cosangle = cose(angle);                                                                                                                                                                
  //  sinangle = sine(angle);                                                                                                                                                                
  cosangle = cos(angle);
  sinangle = sin(angle);
  *re = x*cosangle - y*sinangle;
  *im = x*sinangle + y*cosangle;
}

// N is the number of bins in Pol, which can be a fraction or more than NFreq, which is the number of bins to be correlated.
void CopytoCorrelator(int NFreqBand, int NCorr, float LowerFreqCorr, float UpperFreqCorr, float LowerFreqBand,
                      fftwf_complex *Pol, fftwf_complex *Correlator)
{
  int i, CorrPos;
  for (i=0; i<NFreqBand; i++) {
    CorrPos = (int) (NCorr * (LowerFreqBand-LowerFreqCorr) / (UpperFreqCorr-LowerFreqCorr) + i);
    if (CorrPos>=0 && CorrPos<NCorr) {
      Correlator[CorrPos][0] = Pol[i][0];
      Correlator[CorrPos][1] = Pol[i][1];
    }
  }
}

// Data contains X and Y complex sampled voltages
void Convert2Circular(int N, fftwf_complex *Pol1, fftwf_complex *Pol2, int SwapPols)
{
  int i;
  fftwf_complex *Out1, *Out2;
  Out1 = (fftwf_complex*) fftwf_malloc(N*sizeof(fftwf_complex));
  Out2 = (fftwf_complex*) fftwf_malloc(N*sizeof(fftwf_complex));
#define SQRT22 0.707106781186547
  for (i=0; i<N; i++){
    Out1[i][0] = ( Pol1[i][0]                           + Pol2[i][1] ) * SQRT22;
    Out1[i][1] = (              Pol1[i][1] - Pol2[i][0] ) * SQRT22;
    Out2[i][0] = ( Pol1[i][0]                           - Pol2[i][1] ) * SQRT22;
    Out2[i][1] = (              Pol1[i][1] + Pol2[i][0]  ) *SQRT22;
  }
#undef SQRT22
  for (i=0; i<N; i++){
    Pol1[i][0] = Out1[i][0];
    Pol1[i][1] = Out1[i][1];
    Pol2[i][0] = Out2[i][0];
    Pol2[i][1] = Out2[i][1];
  }
  fftwf_free(Out1);
  fftwf_free(Out2);
}

#define REAL_MLT(R1, I1, R2, I2) (R1)*(R2)-(I1)*(I2)
#define IM_MLT(R1, I1, R2, I2) (R1)*(I2)+(I1)*(R2)

void Calibrate(int NFreq, fftwf_complex *Pol1, fftwf_complex *Pol2, PolCal_type *PolCal)
{
  int i;
  fftwf_complex *Out1, *Out2;
  Out1 = (fftwf_complex*) fftwf_malloc(NFreq*sizeof(fftwf_complex));
  Out2 = (fftwf_complex*) fftwf_malloc(NFreq*sizeof(fftwf_complex));

  for (i=0; i<NFreq; i++){
    Out1[i][0] = REAL_MLT(PolCal->PolCal_11_re[i], PolCal->PolCal_11_im[i], Pol1[i][0], Pol1[i][1]) +
      REAL_MLT(PolCal->PolCal_12_re[i], PolCal->PolCal_12_im[i], Pol2[i][0], Pol2[i][1]);
    Out1[i][1] = IM_MLT(PolCal->PolCal_11_re[i], PolCal->PolCal_11_im[i], Pol1[i][0], Pol1[i][1]) +
      IM_MLT(PolCal->PolCal_12_re[i], PolCal->PolCal_12_im[i], Pol2[i][0], Pol2[i][1]);
    Out2[i][0] = REAL_MLT(PolCal->PolCal_21_re[i], PolCal->PolCal_21_im[i], Pol1[i][0], Pol1[i][1]) +
      REAL_MLT(PolCal->PolCal_22_re[i], PolCal->PolCal_22_im[i], Pol2[i][0], Pol2[i][1]);
    Out2[i][1] = IM_MLT(PolCal->PolCal_21_re[i], PolCal->PolCal_21_im[i], Pol1[i][0], Pol1[i][1]) +
      IM_MLT(PolCal->PolCal_22_re[i], PolCal->PolCal_22_im[i], Pol2[i][0], Pol2[i][1]);
  }
  for (i=0; i<NFreq; i++){
    Pol1[i][0] = Out1[i][0];
    Pol1[i][1] = Out1[i][1];
    Pol2[i][0] = Out2[i][0];
    Pol2[i][1] = Out2[i][1];
  }
  fftwf_free(Out1);
  fftwf_free(Out2);
}
#undef REAL_MLT
#undef IM_MLT

void TransformComplexData(fftwf_complex *FT1, int N, int bandsign)
{
  int i;
  fftwf_complex *FT2;
  FT2 = (fftwf_complex*) fftwf_malloc((N+2)*sizeof(fftwf_complex));
  // Unpack the packed analytic signal
  for (i=0;i<N/2; i++) {
    FT2[i][0] = FT1[i+N/2][0];
    FT2[i][1] = FT1[i+N/2][1];
  }
  for (i=N/2;i<N; i++) {
    FT2[i][0] = FT1[i-N/2][0];
    FT2[i][1] = FT1[i-N/2][1];
  }

  if (bandsign>0)
    for (i=0; i<N; i++) {
      FT1[i][0] = FT2[i][0];
      FT1[i][1] = FT2[i][1];
    }
  else {
    for (i=0; i<N; i++) { // If negative band, invert the band                                                                                                                               
      FT1[i][0] = FT2[N-i][0];
      FT1[i][1] = -FT2[N-i][1]; // Take complex conjugate because band is swapped                                                                                                            
    }
  }
  free(FT2);
}

double CalcRMS(int N, double *Data)
{
  double mean=0, rms=0;
  int i;
  for (i=0; i<N; i++) mean += Data[i]/N;
  for (i=0; i<N; i++) rms += SQR(Data[i]-mean)/N;
  rms=sqrt(rms);
  return rms;
}

// This reads the phase corrections due to the bandpass
void ReadPhaseCorrections(PCF_type *PCF, Fileinfo_type *Fileinfo)
{
  FILE* IO;
  int i, j, Freqchan;
  IO = Sopen(PCF->PhaseCorrectionsFile, "r");
  fscanf(IO, "%d\n", &PCF->PCBins);
  PCF->PhaseCorrections_Pol1 = Make2DArray_float(Fileinfo->NTelescopes-1, PCF->PCBins);
  PCF->PhaseCorrections_Pol2 = Make2DArray_float(Fileinfo->NTelescopes-1, PCF->PCBins);
  for (i=0; i<Fileinfo->NTelescopes-1; i++)
    for (j=0; j<PCF->PCBins; j++) {
      fscanf(IO, "%d %f %f", &Freqchan, &PCF->PhaseCorrections_Pol1[i][j], &PCF->PhaseCorrections_Pol2[i][j]);
      //      PCF->PhaseCorrections_Pol2[i][j]=PCF->PhaseCorrections_Pol1[i][j]; // Uncomment for only 1-pol PCF                                                                             
      if (Freqchan !=j ) { fprintf(stderr, "Error reading %s, line %d\n", PCF->PhaseCorrectionsFile, PCF->PCBins*i+j); exit(0); }
    }
  fclose(IO);
  return;
}

// This routine finds the closest frequency channel given all the Jones matrices and a frequency. It then fills the PolCal with the proper Jones matrix.                                     
// It assumes the list of Jones matrices is incremental in frequency!                                                                                                                        
void FindFreqChannel(int NJonesFreq, int PolFreq, PolCal_type *PolCal, float **Jones, float Freq)
{
  int i, FreqChan;

  if (Freq < Jones[0][0]) FreqChan=0;
  else if (Freq > Jones[NJonesFreq-1][0]) FreqChan=NJonesFreq-1;
  else {
    FreqChan=0;
    for (i=0; i<NJonesFreq-1; i++)
      if (fabs(Freq-Jones[i][0]) > fabs(Freq-Jones[i+1][0])) FreqChan=i+1;
  }
  PolCal->PolCal_11_re[PolFreq] = Jones[FreqChan][1];
  PolCal->PolCal_11_im[PolFreq] = Jones[FreqChan][2];
  PolCal->PolCal_12_re[PolFreq] = Jones[FreqChan][3];
  PolCal->PolCal_12_im[PolFreq] = Jones[FreqChan][4];
  PolCal->PolCal_21_re[PolFreq] = Jones[FreqChan][5];
  PolCal->PolCal_21_im[PolFreq] = Jones[FreqChan][6];
  PolCal->PolCal_22_re[PolFreq] = Jones[FreqChan][7];
  PolCal->PolCal_22_im[PolFreq] = Jones[FreqChan][8];
}

// Read the polarisation calibration for every telescope                                                                                                                                     
void ReadPolCal(PolCal_type *PolCal, Fileinfo_type *Fileinfo)
{
  FILE *IO;
  int i, j;
  float dummy;
  float **Jones;
  int NJonesFreq;
  for (i=0; i<Fileinfo->NTelescopes; i++) {
    if (!PolCal[i].OFF) {
      IO = Sopen(PolCal[i].PolCalFile, "r"); // Open polcalfile and find number of channels                                                                                                  
      NJonesFreq=0;
      while (fscanf(IO, "%f %f %f %f %f %f %f %f %f", &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy) == 9) NJonesFreq++;
      fprintf(stderr, "number of calbins[%d] = %d\n", i, NJonesFreq);
      rewind(IO);
      Jones = Make2DArray_float(NJonesFreq, 9);
      PolCal[i].PolCal_11_re = (float*) calloc(Fileinfo->ReadBins[i], sizeof(float));
      PolCal[i].PolCal_11_im = (float*) calloc(Fileinfo->ReadBins[i], sizeof(float));
      PolCal[i].PolCal_12_re = (float*) calloc(Fileinfo->ReadBins[i], sizeof(float));
      PolCal[i].PolCal_12_im = (float*) calloc(Fileinfo->ReadBins[i], sizeof(float));
      PolCal[i].PolCal_21_re = (float*) calloc(Fileinfo->ReadBins[i], sizeof(float));
      PolCal[i].PolCal_21_im = (float*) calloc(Fileinfo->ReadBins[i], sizeof(float));
      PolCal[i].PolCal_22_re = (float*) calloc(Fileinfo->ReadBins[i], sizeof(float));
      PolCal[i].PolCal_22_im = (float*) calloc(Fileinfo->ReadBins[i], sizeof(float));
      for (j=0; j<NJonesFreq; j++) { // Read polcalfile
        fscanf(IO, "%f %f %f %f %f %f %f %f %f", &Jones[j][0], &Jones[j][1], &Jones[j][2], &Jones[j][3], &Jones[j][4], &Jones[j][5], &Jones[j][6], &Jones[j][7], &Jones[j][8]);
      }
      fclose(IO);
      for (j=0; j<Fileinfo->ReadBins[i]; j++) {
        FindFreqChannel(NJonesFreq, j, &(PolCal[i]), Jones, Fileinfo->Freqinfo[i].LowerFreqBand[0]+Fileinfo->Freqinfo[i].BW * j/Fileinfo->ReadBins[i]);
      }
      Free2DArray_float(Jones, NJonesFreq);
    }
  }
  return;
}

// Read Calc file                                                                                                                                                                            
// The first time this is called, set SetCountCalcSteps to 1. This will then determine how far into the calc file to read.                                                                   
// The next times this routine is called, CountCalcSteps should be increased by 1 before calling and SetCountCalcSteps should be 0                                                           
void ReadCalcfile(Calc_type *Calc, Fileinfo_type *Fileinfo, int SetCountCalcSteps)
{
  int i, j, MJDint, MJDsec;
  double Refdelay, *RefTelPoly;
  FILE *IO;
  char line[MaxChars], Dummy[MaxChars];
  int RefTel = Fileinfo->RefTel;

  RefTelPoly = (double*) calloc(NPolys, sizeof(double));
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
  for (i=0; i<(2*6*Fileinfo->NTelescopes+2)*Calc->CountCalcSteps; i++) fgets(line, MaxChars, IO); // Skip to right timesection                                                               
  for (i=0; i<Fileinfo->NTelescopes; i++) {
    fgets(line, MaxChars, IO);
    sscanf(line, "%23c %lf %lf %lf %lf %lf %lf", Dummy, &Calc->Poly[i][0], &Calc->Poly[i][1], &Calc->Poly[i][2], &Calc->Poly[i][3],
           &Calc->Poly[i][4], &Calc->Poly[i][5]);
    for (j=0; j<5; j++) fgets(line, MaxChars, IO);
  }
  fclose(IO);
  // And again for the reference telescope                                                                                                                                                   
  IO = Sopen(Calc->filename_ref, "r");
  for (i=0; i<22; i++) fgets(line, MaxChars, IO); // Skip first lines                                                                                                                        
  for (i=0; i<(2*6*2+2)*Calc->CountCalcSteps; i++) fgets(line, MaxChars, IO); // Skip to right timesection                                                                                   
  fgets(line, MaxChars, IO);
  sscanf(line, "%23c %lf %lf %lf %lf %lf %lf", Dummy, &RefTelPoly[0], &RefTelPoly[1], &RefTelPoly[2], &RefTelPoly[3], &RefTelPoly[4], &RefTelPoly[5]);
  for (j=0; j<5; j++) fgets(line, MaxChars, IO);
  fclose(IO);

  for (i=0; i<Fileinfo->NTelescopes; i++)
    for (j=0; j<NPolys; j++)
      Calc->Poly[i][j] -= RefTelPoly[j];
  free(RefTelPoly);
}

/* This routine shifts the data <Shift> bins to the right, cyclic. */
void ShiftData(int NBins, double complex **Pol, int FreqChan, int Shift)
{
  int i;
  double complex *Temp;
  Shift = Shift % NBins;
  if (Shift < 0) fprintf(stderr, "ERROR, NEGATIVE SHIFT!\n");
  Temp = (double complex*) calloc(NBins, sizeof(double complex));
  for (i=0; i<NBins; i++) Temp[i] = Pol[i][FreqChan];
  for (i=0; i<Shift; i++) Pol[NBins-Shift+i][FreqChan] = Temp[i];
  for (i=Shift; i<NBins; i++) Pol[i-Shift][FreqChan] = Temp[i];
  free(Temp);
}

/* This function dedisperses the 2-D data */
/* Freq refers to the bottom frequency.*/
void Dedisperse(int NFreqBins, int NTimeBins, double complex **Pol, float FreqMax, float FreqMin, double DM, double BinSize, int PhaseShift)
{
  int i, j, Shift;
  double *FreqValues;
  float Channelwidth = (FreqMax-FreqMin)/NFreqBins;
  FreqValues = (double*) calloc(NFreqBins, sizeof(double));
  for (i=0; i<NFreqBins; i++) FreqValues[i] = FreqMin + (i+0.5)*Channelwidth;
  for (i=0; i<NFreqBins; i++) {
    // Shift = (int) (4149*DM/SQR(FreqValues[i]/1000)/1.E6 / BinSize); /* Shift in bins given by frequency and DM */                                                                         
    Shift = (int) (DM/2.41e-4/SQR(FreqValues[i]/1000.)/1.E6 / BinSize + 0.5) + PhaseShift; /* Shift in bins given by frequency and DM */
    // fprintf(stderr,"%f\t%d\t%12.16f\t%f\n", FreqValues[i], Shift, BinSize, DM);                                                                                                           
    ShiftData(NTimeBins, Pol, i, Shift);                          /* Shift data to correct for (or introduce) dispersion*/
  }
  free(FreqValues);
}

void OutputVis(Fileinfo_type *Fileinfo, double complex ***Vis_pol1, double complex ***Vis_pol2, double complex ***Vis_pol3,
               double complex ***Vis_pol4, int correlate, int NSteps, int NBins, long *DivideBin, double PSRphase)
{
  int i, j, k, l, baseline;
  int Shift;
  if (correlate==1) { // print normal visibilities                                                                                                                                           
    printf("# Timedump Tel1 Tel2 FreqChan Phase1 Amp1 Phase2 Amp2:\n");
    for (i=0; i<(int)((NSteps-DELTA)/Fileinfo->IntegrateNBins)+1; i++) {
      baseline = 0;
      for (j=0; j<Fileinfo->NTelescopes; j++)
        for (k=j; k<Fileinfo->NTelescopes; k++){
          for (l=0; l<Fileinfo->NFreqCorr; l++) {
            if (DOAC || j!=k)
              if (l==0) // Put the very edge of the band to the value of its neighbour                                                                                                       
                if (COMPLEXVIS)
                  printf("%d %d %d %d %f %f %f %f %f %f %f %f\n", i, j, k, l, creal(Vis_pol1[i][baseline][1]), cimag(Vis_pol1[i][baseline][1]),
                         creal(Vis_pol2[i][baseline][1]), cimag(Vis_pol2[i][baseline][1]), creal(Vis_pol3[i][baseline][1]), cimag(Vis_pol3[i][baseline][1]),
                         creal(Vis_pol4[i][baseline][1]), cimag(Vis_pol4[i][baseline][1]));
                else
                  printf("%d %d %d %d %f %f %f %f %f %f %f %f\n", i, j, k, l, carg(Vis_pol1[i][baseline][1]), cabs(Vis_pol1[i][baseline][1]),
                         carg(Vis_pol2[i][baseline][1]), cabs(Vis_pol2[i][baseline][1]), carg(Vis_pol3[i][baseline][1]), cabs(Vis_pol3[i][baseline][1]),
                         carg(Vis_pol4[i][baseline][1]), cabs(Vis_pol4[i][baseline][1]));
              else
                if (COMPLEXVIS)
                  printf("%d %d %d %d %f %f %f %f %f %f %f %f\n", i, j, k, l, creal(Vis_pol1[i][baseline][l]), cimag(Vis_pol1[i][baseline][l]),
                         creal(Vis_pol2[i][baseline][l]), cimag(Vis_pol2[i][baseline][l]), creal(Vis_pol3[i][baseline][l]), cimag(Vis_pol3[i][baseline][l]),
                         creal(Vis_pol4[i][baseline][l]), cimag(Vis_pol4[i][baseline][l]));
                else
                  printf("%d %d %d %d %f %f %f %f %f %f %f %f\n", i, j, k, l, carg(Vis_pol1[i][baseline][l]), cabs(Vis_pol1[i][baseline][l]),
                         carg(Vis_pol2[i][baseline][l]), cabs(Vis_pol2[i][baseline][l]), carg(Vis_pol3[i][baseline][l]), cabs(Vis_pol3[i][baseline][l]),
                         carg(Vis_pol4[i][baseline][l]), cabs(Vis_pol4[i][baseline][l]));
          }
          if (DOAC || j!=k) {
            printf("\n");
            baseline++;
          }
        }
    }
  }
  if (correlate==2) { // print pulsar profile                                                                                                                                                
    for (k=0; k<Fileinfo->NBaselines; k++)
      for (i=0; i<NBins; i++) {
        Vis_pol1[k][i][0] = Vis_pol1[k][i][1];
        Vis_pol2[k][i][0] = Vis_pol2[k][i][1];
        Vis_pol3[k][i][0] = Vis_pol3[k][i][1];
        Vis_pol4[k][i][0] = Vis_pol4[k][i][1];
        for (j=0; j<Fileinfo->NFreqCorr; j++) {
          Vis_pol1[k][i][j] /= DivideBin[i];
          Vis_pol2[k][i][j] /= DivideBin[i];
          Vis_pol3[k][i][j] /= DivideBin[i];
          Vis_pol4[k][i][j] /= DivideBin[i];
        }
      }
    // Dedisperse every baseline                                                                                                                                                             
    baseline=0;
    Shift = (int) NBins - (PSRphase*NBins) + NBins/2;
    printf("# Shifting %d bins\n", Shift);
    if (Shift < 0) fprintf(stderr, "ERROR, NEGATIVE SHIFT in outputvis!\n");
    for (j=0; j<Fileinfo->NTelescopes; j++)
      for (k=j; k<Fileinfo->NTelescopes; k++)
        if (DOAC || j!=k) {
          Dedisperse(Fileinfo->NFreqCorr, NBins, Vis_pol1[baseline], Fileinfo->UpperFreqCorr, Fileinfo->LowerFreqCorr, Fileinfo->DM, Fileinfo->IntegrateSeconds, Shift);
          Dedisperse(Fileinfo->NFreqCorr, NBins, Vis_pol2[baseline], Fileinfo->UpperFreqCorr, Fileinfo->LowerFreqCorr, Fileinfo->DM, Fileinfo->IntegrateSeconds, Shift);
          Dedisperse(Fileinfo->NFreqCorr, NBins, Vis_pol3[baseline], Fileinfo->UpperFreqCorr, Fileinfo->LowerFreqCorr, Fileinfo->DM, Fileinfo->IntegrateSeconds, Shift);
          Dedisperse(Fileinfo->NFreqCorr, NBins, Vis_pol4[baseline], Fileinfo->UpperFreqCorr, Fileinfo->LowerFreqCorr, Fileinfo->DM, Fileinfo->IntegrateSeconds, Shift);
          baseline++;
        }
    for (i=0; i<NBins+1; i++){
      baseline = 0;
      for (j=0; j<Fileinfo->NTelescopes; j++)
        for (k=j; k<Fileinfo->NTelescopes; k++){
          for (l=0; l<Fileinfo->NFreqCorr; l++)
            if (DOAC || j!=k) {
              printf("%d %d %d %d %f %f %f %f %f %f %f %f\n", i, j, k, l, creal(Vis_pol1[baseline][i][l]), cimag(Vis_pol1[baseline][i][l]),
                     creal(Vis_pol2[baseline][i][l]), cimag(Vis_pol2[baseline][i][l]), creal(Vis_pol3[baseline][i][l]), cimag(Vis_pol3[baseline][i][l]),
                     creal(Vis_pol4[baseline][i][l]), cimag(Vis_pol4[baseline][i][l]));
            }
          if (DOAC || j!=k) {
            printf("\n");
            baseline++;
          }
        }
    }
  }
  if (correlate==3) { // print time-series                                                                                                                                                   
    baseline = 0;
    for (j=0; j<Fileinfo->NTelescopes; j++)
      for (k=j; k<Fileinfo->NTelescopes; k++)
        if (DOAC || j!=k) {
          Dedisperse(Fileinfo->NFreqCorr, NSteps/Fileinfo->IntegrateNBins, Vis_pol1[baseline], Fileinfo->UpperFreqCorr, Fileinfo->LowerFreqCorr, Fileinfo->DM,
                     Fileinfo->IntegrateSeconds, PSRphase);
          Dedisperse(Fileinfo->NFreqCorr, NSteps/Fileinfo->IntegrateNBins, Vis_pol2[baseline], Fileinfo->UpperFreqCorr, Fileinfo->LowerFreqCorr, Fileinfo->DM,
                     Fileinfo->IntegrateSeconds, PSRphase);
          baseline++;
        }
    for (i=0; i<(int)((NSteps-DELTA)/Fileinfo->IntegrateNBins)+1; i++) {
      baseline = 0;
      for (j=0; j<Fileinfo->NTelescopes; j++)
        for (k=j; k<Fileinfo->NTelescopes; k++) {
          for (l=0; l<Fileinfo->NFreqCorr; l++)
            if (DOAC || j!=k) printf("%d %d %d %d %f %f %f %f\n", i, j, k, l, carg(Vis_pol1[baseline][i][l]), cabs(Vis_pol1[baseline][i][l]),
                                     carg(Vis_pol2[baseline][i][l]), cabs(Vis_pol2[baseline][i][l]));
          if (DOAC || j!=k) {
            baseline++;
            printf("\n");
          }
        }
    }
  }
}


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
          //      Fileinfo->ClockDrift[i] += (ClockOffset_read - ClockOffset_read_prev)/86400*(Fileinfo->LowerFreqCorr+0.5*Fileinfo->BWCorr)*2*PI; // in rad/s                               
          //      printf("# Clockdrift[%d] = %f rad/s\n", i, Fileinfo->ClockDrift[i]);                                                                                                       
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
// This routine gets the geometrical offset of the start of the correlation. This is applied as soon as the files are opened.                                                                
double GetStartGeoOffset(Calc_type Calc, Fileinfo_type *Fileinfo, int Telescope, long reducedbinnumber)
{
  int i;
  long binoffset;
  double timedelay = 0;
  for (i=0; i<NPolys; i++) {
    timedelay += pow(reducedbinnumber*Fileinfo->tsamp, i) * Calc.Poly[Telescope][i]; // in us. Gives timedelay from 0 to now                                                                 
  }
  return -timedelay*1000; // convert to ns                                                                                                                                                   
}

// Calculate the full geometric timedelay of one Telescope in units of us                                                                                                                    
// reducedbinnumber is the binnumber from the start of the current calc-step to the current bin in units of output bins                                                                      
double GetTimedelay(int N, int Telescope, Calc_type *Calc, Fileinfo_type *Fileinfo, long reducedbinnumber)
{
  int i;
  double timedelay = 0;
  for (i=0; i<NPolys; i++)
    timedelay += pow(reducedbinnumber*Fileinfo->tsamp, i) * Calc->Poly[Telescope][i]; // in us. Gives timedelay from 0 to now                                                                
  timedelay += Fileinfo->GeoShifted[Telescope]*Fileinfo->SampTime[Telescope]/1000; // subtract bins that are already shifted                                                                 
  return timedelay; // in us                                                                                                                                                                 
}


// Calculate the fringerotation                                                                                                                                                              
// time is the time from the start of the calcfile to now in seconds                                                                                                                         
// timedelay is in us                                                                                                                                                                        
double GetFringeRotation(int N, int Telescope, int freqcounter, double timedelay, Calc_type *Calc, Fileinfo_type *Fileinfo, double time, int Pol)
{
  double FringeRotation, ReducedFringeRotation, fracdelay, freqoffset;
  if (Pol==0) {
    fracdelay = Fileinfo->Fracdelay_pol1[Telescope];   // in us                                                                                                                              
    freqoffset = Fileinfo->Freqoffset_pol1[Telescope]; // in radians                                                                                                                         
  }
  else {
    fracdelay = Fileinfo->Fracdelay_pol2[Telescope];   // in us                                                                                                                              
    freqoffset = Fileinfo->Freqoffset_pol2[Telescope]; // in radians                                                                                                                         
  }

  // Correct for Fringerotation due to binshifting the clockoffset and geometric delay                                                                                                       
  if (Fileinfo->IsComplex[Telescope])
    FringeRotation = -2*PI*(Fileinfo->BinsShifted[Telescope] - (timedelay+fracdelay)*1000/Fileinfo->SampTime[Telescope])
      * (Fileinfo->Freqinfo[Telescope].LowerFreqBand[0] + Fileinfo->Freqinfo[Telescope].BW/2) / Fileinfo->Freqinfo[Telescope].BW;
  else
    FringeRotation = -2*PI*(Fileinfo->BinsShifted[Telescope] - (timedelay+fracdelay)*1000/Fileinfo->SampTime[Telescope])
      * (Fileinfo->Freqinfo[Telescope].LowerFreqBand[0]) / Fileinfo->Freqinfo[Telescope].BW;

  // Add correction for the Fractional Sample Error Correction + atmospheric delay                                                                                                           
  if (Fileinfo->IsComplex[Telescope])
    FringeRotation += 2*PI*(timedelay+fracdelay)*((freqcounter-N/2)/(float)N*Fileinfo->BWCorr) + freqoffset - time * Fileinfo->FringeDrift[Telescope] - time * Fileinfo->ClockDrift[Telescop\
e] *
      ((freqcounter-N/2)/(float)N*Fileinfo->BWCorr + Fileinfo->Freqinfo[Telescope].LowerFreqBand[0]) * 2*PI;
  else
    FringeRotation += 2*PI*(timedelay+fracdelay)*(freqcounter/(float)N*Fileinfo->BWCorr) + freqoffset - time * Fileinfo->FringeDrift[Telescope] - time * Fileinfo->ClockDrift[Telescope] *
      (freqcounter/(float)N*Fileinfo->BWCorr + Fileinfo->Freqinfo[Telescope].LowerFreqBand[0]) * 2*PI;

  return fmod(FringeRotation, 2*PI);
}


// This routine aligns the data from every telescope before addition or multiplication. The reference is always telescope 0.                                                                 
// Think about reference telescope !=0 !                                                                                                                                                     
void AlignData(int N, PCF_type *PCF, FFT_type *FFT, Calc_type *Calc, Fileinfo_type *Fileinfo, long reducedbinnumber,
               int RefreshRotation, double time, int countloop)
{
  int i, j;
  double timedelay;
  for (j=1; j<Fileinfo->NTelescopes; j++) { // Loop over telescopes, skip telescope 0                                                                                                        
    if (RefreshRotation && j) timedelay = GetTimedelay(N, j, Calc, Fileinfo, reducedbinnumber);  //timedelay in us                                                                           
    if (Fileinfo->IncoherentAddition == 1) { // Add incoherently (meaning: add an extra phaseshift that is constant over frequency and varies linearly over time)                            
      for (i=1; i<N; i++) { // Loop over frequency channels                                                                                                                                  
        if (RefreshRotation) {
          Fileinfo->FringeRotationPol1[j][i] =
            PCF->PhaseCorrections_Pol1[j-1][(int)(i/(float)N*PCF->PCBins)] - GetFringeRotation(N, j, i, timedelay, Calc, Fileinfo, time, 0);
          Fileinfo->FringeRotationPol2[j][i] =
            PCF->PhaseCorrections_Pol2[j-1][(int)(i/(float)N*PCF->PCBins)] - GetFringeRotation(N, j, i, timedelay, Calc, Fileinfo, time, 1);
        }
        // Align vectors with telescope 0 and add incoherent phase                                                                                                                           
        RotateVec(&FFT->Correlate_pol1[countloop][j][i][0], &FFT->Correlate_pol1[countloop][j][i][1], Fileinfo->FringeRotationPol1[j][i] + j*Fileinfo->IncoherentAddition_Phase);
        RotateVec(&FFT->Correlate_pol2[countloop][j][i][0], &FFT->Correlate_pol2[countloop][j][i][1], Fileinfo->FringeRotationPol2[j][i] + j*Fileinfo->IncoherentAddition_Phase);
      }
      Fileinfo->IncoherentAddition_Phase += IncoherentPhaseStep; // Increase the phase for the next time-sample                                                                              
      //      fprintf(stderr, "IncoherentAddition_Phase = %f\n", Fileinfo->IncoherentAddition_Phase);                                                                                        
    }
    else // Add coherently                                                                                                                                                                   
      for (i=1; i<N; i++) { // Loop over frequency channels, don't rotate bin 0, as it only contains the DC peak and should therefore be real
        if (RefreshRotation) {
          Fileinfo->FringeRotationPol1[j][i] =
            PCF->PhaseCorrections_Pol1[j-1][(int)(i/(float)N*PCF->PCBins)] - GetFringeRotation(N, j, i, timedelay, Calc, Fileinfo, time, 0);
          Fileinfo->FringeRotationPol2[j][i] =
            PCF->PhaseCorrections_Pol2[j-1][(int)(i/(float)N*PCF->PCBins)] - GetFringeRotation(N, j, i, timedelay, Calc, Fileinfo, time, 1);
        }
        // Align vectors with telescope 0                                                                                                                                                    
        RotateVec(&FFT->Correlate_pol1[countloop][j][i][0], &FFT->Correlate_pol1[countloop][j][i][1], Fileinfo->FringeRotationPol1[j][i]);
        RotateVec(&FFT->Correlate_pol2[countloop][j][i][0], &FFT->Correlate_pol2[countloop][j][i][1], Fileinfo->FringeRotationPol2[j][i]);
      }
  } // End loop over telescopes
}

void AddtoSum(int N, fftwf_complex *sum_pol1, fftwf_complex *sum_pol2, FFT_type *FFT, Fileinfo_type *Fileinfo, Bandpass_type *Bandpass, Zap_type *Zap, int countloop)
{
  int i, j;
  float divider;

  for (i=0; i<N; i++) { // Loop over frequency channels                                                                                                                                      
    sum_pol1[i][0] = 0;
    sum_pol1[i][1] = 0;
    sum_pol2[i][0] = 0;
    sum_pol2[i][1] = 0;
  }

  for (j=0; j<Fileinfo->NTelescopes; j++) { // Loop over telescopes                                                                                                                          
    if (j == SUMONETELESCOPE || SUMONETELESCOPE == -1) {
      for (i=0; i<N; i++) { // Loop over frequency channels                                                                                                                                  
        if (Zap->Zap[j][MatchInt(N,i,Zap->ZapBins)]) {
          divider=Bandpass->Bandpass[0][MatchInt(N,i,Bandpass->BandpassBins)]/Bandpass->Bandpass[j][MatchInt(N,i,Bandpass->BandpassBins)];
          // IF ZERO-DM IS ON, WE SHOULD NOT DIVIDE BY RMS                                                                                                                                   
          sum_pol1[i][0] += Fileinfo->Scale_Pol1[j]*FFT->Correlate_pol1[countloop][j][i][0]/Fileinfo->rms[j].pol1 / divider;
          sum_pol1[i][1] += Fileinfo->Scale_Pol1[j]*FFT->Correlate_pol1[countloop][j][i][1]/Fileinfo->rms[j].pol1 / divider;
          sum_pol2[i][0] += Fileinfo->Scale_Pol2[j]*FFT->Correlate_pol2[countloop][j][i][0]/Fileinfo->rms[j].pol2 / divider;
          sum_pol2[i][1] += Fileinfo->Scale_Pol2[j]*FFT->Correlate_pol2[countloop][j][i][1]/Fileinfo->rms[j].pol2 / divider;
        }
      }
    }
  }
}

// reducedbinnumber is in units of 50 * NFREQ ns                                                                                                                                             
// Fileinfo->Fracdelay is in units of ns, but gets converted to us                                                                                                                           
void AddtoVis(int N, int Telescope1, int Telescope2, double complex *vis_pol1, double complex *vis_pol2, double complex *vis_pol3, double complex *vis_pol4, FFT_type *FFT,
              Fileinfo_type *Fileinfo, int countloop)
{
  int i;
  complex p1t10,p1t11,p1t20,p1t21,p2t10,p2t11,p2t20,p2t21;
  if (Telescope1 != Telescope2) {
    //    for (i=0; i<N; i++) {                                                                                                                                                              
    for (i=1; i<N; i++) {
      p1t10 = FFT->Correlate_pol1[countloop][Telescope1][i][0];
      p1t11 = FFT->Correlate_pol1[countloop][Telescope1][i][1];
      p1t20 = FFT->Correlate_pol1[countloop][Telescope2][i][0];
      p1t21 = FFT->Correlate_pol1[countloop][Telescope2][i][1];
      p2t10 = FFT->Correlate_pol2[countloop][Telescope1][i][0];
      p2t11 = FFT->Correlate_pol2[countloop][Telescope1][i][1];
      p2t20 = FFT->Correlate_pol2[countloop][Telescope2][i][0];
      p2t21 = FFT->Correlate_pol2[countloop][Telescope2][i][1];
      // Calculate polarisation                                                                                                                                                              
      vis_pol1[i] += (p1t10*p1t20 + p1t11*p1t21) + I * (p1t11*p1t20 - p1t10*p1t21);
      vis_pol2[i] += (p2t10*p2t20 + p2t11*p2t21) + I * (p2t11*p2t20 - p2t10*p2t21);
      // And cross polarisation                                                                                                                                                              
      vis_pol3[i] += (p1t10*p2t20 + p1t11*p2t21) + I * (p1t11*p2t20 - p1t10*p2t21);
      vis_pol4[i] += (p2t10*p1t20 + p2t11*p1t21) + I * (p2t11*p1t20 - p2t10*p1t21);
    }
  }
  else // Do Auto Correlation                                                                                                                                                                
    for (i=1; i<N; i++) {
      p1t10 = FFT->Correlate_pol1[countloop][Telescope1][i][0];
      p1t11 = FFT->Correlate_pol1[countloop][Telescope1][i][1];
      p1t20 = FFT->Correlate_pol1[countloop][Telescope2][i][0];
      p1t21 = FFT->Correlate_pol1[countloop][Telescope2][i][1];
      p2t10 = FFT->Correlate_pol2[countloop][Telescope1][i][0];
      p2t11 = FFT->Correlate_pol2[countloop][Telescope1][i][1];
      p2t20 = FFT->Correlate_pol2[countloop][Telescope2][i][0];
      p2t21 = FFT->Correlate_pol2[countloop][Telescope2][i][1];
      vis_pol1[i] += Power(p1t10, p1t11); // LL                                                                                                                                              
      vis_pol2[i] += Power(p2t10, p2t11); // RR                                                                                                                                              
      vis_pol3[i] += (p1t10*p2t10 + p1t11*p2t11) + I * (p1t11*p2t10 - p1t10*p2t11);
      vis_pol4[i] += (p2t10*p1t10 + p2t11*p1t11) + I * (p2t11*p1t10 - p2t10*p1t11);
    }
}

