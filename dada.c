#include "timeseries.h"

#define HEADERSIZE 4096

// Read an unsigned 64bit integer value from a character buffer
uint64_t read_uint64_value(char *keyword,char *buffer)
{
  uint64_t result;
  char *pch;
  
  if ((pch=strstr(buffer,keyword))!=NULL)
    sscanf(pch+strlen(keyword)+1,"%Lu",&result);
  else {
    fprintf(stderr,"Error reading keyword %s\n",keyword);
    exit;
  }

  return result;
}

// Read an integer value from a character buffer
int read_integer_value(char *keyword,char *buffer)
{
  int result;
  char *pch;

  if ((pch=strstr(buffer,keyword))!=NULL)
    sscanf(pch+strlen(keyword)+1,"%d",&result);
  else {
    fprintf(stderr,"Error reading keyword %s\n",keyword);
    exit;
  }

  return result;
}

// Read a string from a character buffer
char *read_string_value(char *keyword,char *buffer,int len)
{
  char *result;
  char *pch;

  result=(char *) malloc(len*sizeof(char));
  if ((pch=strstr(buffer,keyword))!=NULL) {
    sscanf(pch+strlen(keyword)+1,"%s",result);
  } else {
    fprintf(stderr,"Error reading keyword %s\n",keyword);
    exit;
  }

  return result;
}

// Read an array of doubles from a character buffer
double *read_double_array(char *keyword,char *buffer,int n)
{
  double *result;
  char *pch;

  result=(double *) malloc(sizeof(double)*n);
  if ((pch=strstr(buffer,keyword))!=NULL && n==6) {
    sscanf(pch+strlen(keyword)+1,"%lf %lf %lf %lf %lf %lf",&result[0],&result[1],&result[2],&result[3],&result[4],&result[5]);
  } else {
    fprintf(stderr,"Error reading keyword %s\n",keyword);
    exit;
  }

  return result;
}

// Read a double from a character buffer
double read_double_value(char *keyword,char *buffer)
{
  double result;
  char *pch;

  if ((pch=strstr(buffer,keyword))!=NULL)
    sscanf(pch+strlen(keyword)+1,"%lf",&result);
  else {
    fprintf(stderr,"Error reading keyword %s\n",keyword);
    exit;
  }

  return result;
}

// Read a timestamp
int read_timestamp(char *keyword,char *buffer)
{
  int result;
  int iyear,imonth,iday,ihour,imin,isec;
  char *pch;
  
  if ((pch=strstr(buffer,keyword))!=NULL) {
    sscanf(pch+strlen(keyword)+1,"%04d-%02d-%02d-%02d:%02d:%02d",&iyear,&imonth,&iday,&ihour,&imin,&isec);
  } else {
    fprintf(stderr,"Error reading keyword %s\n",keyword);
    exit;
  }
  result=3600*ihour+60*imin+isec;
  
  return result;
}

// Read DADA header
struct timeseries read_dada_header(FILE *file)
{
  char header[HEADERSIZE];
  struct timeseries s;
  int bytes_read;

  // Read 4096 bytes containing the header
  bytes_read=fread(header,1,HEADERSIZE,file);

  // Get source name
  strcpy(s.source,read_string_value("SOURCE",header,20));
  
  // Get telescope name
  strcpy(s.telescope,read_string_value("TELESCOPE",header,20));

  // Get instrument name
  strcpy(s.instrument,read_string_value("INSTRUMENT",header,20));

  // Get frequency
  s.freq=read_double_value("FREQ",header);

  // Get bandwidth
  s.bw=read_double_value("BW",header);

  // Get sampling time in s
  s.tsamp=read_double_value("TSAMP",header)/1000000.0;

  // Get file_size
  s.file_size=read_uint64_value("FILE_SIZE",header);

  // Get obs_offset
  s.obs_offset=read_uint64_value("OBS_OFFSET",header);

  // Get MJD start
  s.mjd_start=read_double_value("MJD_START",header);
  s.intmjd=(int) s.mjd_start;

  // Get UTC start
  s.intsec=read_timestamp("UTC_START",header);

  // get NBIT, NDIM, NCHAN NPOL
  s.nbit=read_integer_value("NBIT",header);
  s.ndim=read_integer_value("NDIM",header);
  s.npol=read_integer_value("NPOL",header);

  // Some old WSRT files don't have NCHAN specified.
  if (strstr(header,"NCHAN")!=NULL)
    s.nchan=read_integer_value("NCHAN",header);
  else
    s.nchan=1;

  // Get bytes per sample
  s.bytes_per_sample=s.nbit*s.npol*s.nchan*s.ndim/8;

  // Get samples per second
  s.samples_per_second=(uint64_t) ceil(1.0/s.tsamp);

  // Get bytes per second
  s.bytes_per_second=s.samples_per_second*s.bytes_per_sample;

  return s;
}
