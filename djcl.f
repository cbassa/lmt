      SUBROUTINE sla_DJCL (DJM, IY, IM, ID, FD, J)
*+
*     - - - - -
*      D J C L
*     - - - - -
*
*  Modified Julian Date to Gregorian year, month, day,
*  and fraction of a day.
*
*  Given:
*     DJM      dp     modified Julian Date (JD-2400000.5)
*
*  Returned:
*     IY       int    year
*     IM       int    month
*     ID       int    day
*     FD       dp     fraction of day
*     J        int    status:
*                       0 = OK
*                      -1 = unacceptable date (before 4701BC March 1)
*
*  The algorithm is derived from that of Hatcher 1984
*  (QJRAS 25, 53-55).
*
*  P.T.Wallace   Starlink   27 April 1998
*
*  Copyright (C) 1998 Rutherford Appleton Laboratory
*-

      IMPLICIT NONE

      DOUBLE PRECISION DJM
      INTEGER IY,IM,ID
      DOUBLE PRECISION FD
      INTEGER J

      DOUBLE PRECISION F,D
      INTEGER JD,N4,ND10



*  Check if date is acceptable
      IF (DJM.LE.-2395520D0.OR.DJM.GE.1D9) THEN
         J=-1
      ELSE
         J=0

*     Separate day and fraction
         F=MOD(DJM,1D0)
         IF (F.LT.0D0) F=F+1D0
         D=ANINT(DJM-F)

*     Express day in Gregorian calendar
         JD=NINT(D)+2400001

         N4=4*(JD+((6*((4*JD-17918)/146097))/4+1)/2-37)
         ND10=10*(MOD(N4-237,1461)/4)+5

         IY=N4/1461-4712
         IM=MOD(ND10/306+2,12)+1
         ID=MOD(ND10,306)/10+1
         FD=F

         J=0

      END IF

      END
