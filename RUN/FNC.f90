!
!     QSAT:  saturation water vapour specific humidity
!
FUNCTION FQSAT(TT,PP)

  include "./param.inc"

  REAL*8 FQSAT
  REAL*8     TT, PP, QS

!  print *,"IN",PP,TT,TMELT,RVAP
  FQSAT   = EPSV * ES0 / PP  * DEXP(   (EL+EMELT/2.D0*(1.D0-SIGN(1.D0,TT-TQICE)))   /RVAP *( 1.D0/TMELT - 1.D0/TT )           )


END FUNCTION FQSAT


FUNCTION SINB(DOY,HOUR,LAT)

  include "./param.inc"

  REAL*8 SINB
  INTEGER DOY
  REAL*8 HOUR
  REAL*8 LAT

  REAL*8 DEL
  REAL*8 H
  REAL*8 LATRAD


!  DEL = -23.4D0 * PI * COS(2.D0*PI*DBLE(DOY+10)/365.D0) / 180.D0
  DEL = -DASIN(DSIN(23.45D0 * 2.D0 * PI / 360.D0) * DCOS(2.D0 * PI * (DBLE(DOY) + 10.D0)/365.D0))  ! Goudriaan and van Laar (1994)
  H = PI * (HOUR - 12.D0) / 12.D0
  LATRAD = 2.D0 * PI * LAT / 360.D0
  SINB = MAX(0.D0,DSIN(LATRAD) * DSIN(DEL) + DCOS(LATRAD) * DCOS(DEL) * DCOS(H))  ! Solor elevation

END FUNCTION SINB

FUNCTION CALHOUR(IHOUR,TRES)

  REAL*8 CALHOUR
  INTEGER IHOUR
  INTEGER TRES

  CALHOUR =  DBLE(IHOUR-1) * DBLE(TRES) / 60.D0 / 60.D0 + DBLE(TRES)*0.5D0 / 60.D0 / 60.D0

END FUNCTION CALHOUR

FUNCTION DAYL(DOY,LAT)

  include "./param.inc"

  REAL*8 DAYL
  INTEGER DOY
  REAL*8 LAT

  REAL*8 DEL
  REAL*8 LATRAD

  REAL*8 DD

  DEL = -DASIN(DSIN(23.45D0 * 2.D0 * PI / 360.D0) * DCOS(2.D0 * PI * (DBLE(DOY) + 10.D0)/365.D0))  ! Goudriaan and van Laar (1994)
  LATRAD = 2.D0 * PI * LAT / 360.D0

!  print *,DEL,LATRAD
!  print *,DSIN(DEL),DSIN(LATRAD),DCOS(LATRAD),DCOS(DEL)
!  print *,DSIN(DEL)*DSIN(LATRAD)/(DCOS(LATRAD)*DCOS(DEL))
  
  DD=DSIN(LATRAD)*DSIN(DEL) / (DCOS(LATRAD)*DCOS(DEL))

  IF(DD>1.D0)THEN
     DAYL=24.D0
  ELSE IF(DD< -1.D0)THEN
     DAYL=0.D0
  ELSE
     DAYL = 12.D0 * ( 1.D0 + (2.D0/PI)*DASIN( DSIN(LATRAD)*DSIN(DEL) / (DCOS(LATRAD)*DCOS(DEL)) )   )
  END IF

!  print *,DAYL

END FUNCTION DAYL


