SUBROUTINE CROP(YLD,PRM,OCDI,PLT,EMR,GRN,DVS,DVSL,aGDH,aVD,CDI,HDI,TAVE,TMX,NHED,aO3FLX,LAI,LAImx,HGT,ROT,WSH,WSO,WST,WLF,WRT,WAR,WSR,WDL,WGR,SWGR,NGR,NSP,GLF,GST,GRT,GSO,SLN,NFERT,Vsat,GPP,RSP,TMP,PLTDOY1,DOY,HOUR,TRES,GDHm1,hDVS,CFLF,CFST,CFRT,CFSO,VN,TB,TO,TH,LEFY0,LEFX1,LEFY1,LEFX2,LEFY2,LEFX3,LEFY3,PNCLX1,PNCLY1,PNCLX2,PNCLY2,PNCLX3,PNCLY3,DLFX1,DLFY1,DLFX2,DLFY2,DLFX3,DLFY3,RTX,RTY,RTX2,FSTR,SLWYA,SLWYB,SLWX,HGTAA,HGTAB,HGTBA,HGTBB,GZRT,MXRT,SLNX1,SLNY1,SLNX2,SLNX3,SLNYMX,SLNYMN,SLNK,WSTRS,ILON,ILAT,IRR,TCmin,THcrit,HI,TAVE5DAY,TAVE5CNT,HVT_TAVE,LLFst,kLLF,DLFX,CO2PPM,YEAR)

  IMPLICIT NONE

  include "./param.inc"

!  [OUTPUT]
  REAL*8 YLD
  REAL*8 PRM
  REAL*8 OCDI

  REAL*8 LAI    ! LAI [m**2/m**2]
  REAL*8 LAImx
  REAL*8 HGT    ! Height [m]


  REAL*8 SLN
  
  REAL*8 aO3FLX    !OZN




!  [INPUT]
  INTEGER PLTDOY1
  INTEGER DOY
  REAL*8  HOUR



  REAL*8 TMP ![K]
  REAL*8 GPP
  REAL*8 RSP
  INTEGER TRES

  REAL*8 GDHm1
!  REAL*8 GDHm2
  REAL*8 hDVS


  REAL*8  CFLF         ! Fraction of C for leave
  REAL*8  CFST         ! Fraction of C for stem
  REAL*8  CFRT         ! Fraction of C for root
  REAL*8  CFSO         ! Fraction of C for storage orga
  
  INTEGER VN
  REAL*8  TB
  REAL*8  TO
  REAL*8  TH

  REAL*8  RTX
  REAL*8  RTY
  REAL*8  RTX2
  REAL*8  LEFY0
  REAL*8  LEFX1
  REAL*8  LEFY1
  REAL*8  LEFX2
  REAL*8  LEFY2
  REAL*8  LEFX3
  REAL*8  LEFY3
  REAL*8  PNCLX1
  REAL*8  PNCLY1
  REAL*8  PNCLX2
  REAL*8  PNCLY2
  REAL*8  PNCLX3
  REAL*8  PNCLY3
  REAL*8  DLFX1
  REAL*8  DLFY1
  REAL*8  DLFX2
  REAL*8  DLFY2
  REAL*8  DLFX3
  REAL*8  DLFY3

  REAL*8 LLFst
  REAL*8 kLLF
  REAL*8 DLFX


  REAL*8  FSTR
  REAL*8  SLWYA
  REAL*8  SLWYB
  REAL*8  SLWX
  REAL*8  HGTAA
  REAL*8  HGTAB
  REAL*8  HGTBA
  REAL*8  HGTBB
  REAL*8  GZRT        !Growth rate of root [m/day]
  REAL*8  MXRT        !Maximum root length [m]

  REAL*8 SLNX1
  REAL*8 SLNY1
  REAL*8 SLNX2
  REAL*8 SLNX3
  REAL*8 SLNYMX
  REAL*8 SLNYMN
  REAL*8 SLNK

  REAL*8 TCmin
  REAL*8 THcrit

  REAL*8 HI

  REAL*8  WSTRS

  INTEGER ILAT,ILON

  INTEGER IRR

  REAL*8 NFERT
  REAL*8 Vsat

  REAL*8 TAVE5DAY(86400/TRES*5)
  INTEGER TAVE5CNT
  REAL*8 HVT_TAVE

  INTEGER PLANTDOY
  REAL*8 CO2PPM

!  [INTERNAL VARIABLE]



  REAL*8 WSH    ! Biomass living
  REAL*8 WSO    ! Biomass storage organ
  REAL*8 WST    ! Biomass stem
  REAL*8 WLF    ! Biomass leaf
  REAL*8 WRT    ! Biomass root
  REAL*8 WAR    ! Biomass
  REAL*8 WSR    ! Biomass shielded reserve
  REAL*8 WDL    ! Biomass dead leaf
  REAL*8 WGR
  REAL*8 SWGR

  REAL*8 NGR
  REAL*8 NSP

  REAL*8 ROT    ! Root length [m]

  REAL*8 DVS
  REAL*8 DVSL

  REAL*8 ASMN   !net assimilation

  INTEGER PLT
  LOGICAL EMR
  INTEGER GRN

  REAL*8 LLF
  REAL*8 LRT
  REAL*8 LSTR

  REAL*8 aGDH
  REAL*8 aVD

  REAL*8 CDI
  REAL*8 HDI
  REAL*8 TAVE
  REAL*8 TMX
  REAL*8 NHED
  
  REAL*8 GLF
  REAL*8 GST
  REAL*8 GSO
  REAL*8 GRT
  REAL*8 GSR
  REAL*8 GGR
  REAL*8 GWAR

  REAL*8 TAVE5
  INTEGER II

  INTEGER YEAR

  TAVE5 = 0.D0
  DO II=1,TAVE5CNT
     TAVE5 = TAVE5 + TAVE5DAY(II)
  END DO
  TAVE5 = TAVE5 / DBLE(TAVE5CNT)

  PLANTDOY=PLTDOY1
  if(PLANTDOY<1)then
     PLANTDOY = PLTDOY1 + 365
  end if
  if(PLANTDOY>365)then
     PLANTDOY = PLTDOY1 - 365
  end if



  CALL JUDPLT ( PLT,EMR,GRN,DVS,DVSL,WSH,WLF,WST,WSO,WRT,WSR,WAR,WDL,SWGR,WGR,NGR,NSP,LAI,HGT,ROT,aGDH,aO3FLX,CDI,HDI,NHED,DOY,HOUR,PLANTDOY,IRR,DLFX )  


  IF(PLT .GT. 0 )THEN

!     CALL CRODVS( DVS, aGDH,TMP,GDHm1,GDHm2,TRES,PLT  )
     CALL CRODVS( DVS, aGDH,aVD,TMP,GDHm1,TRES,PLT,VN,TB,TO,TH ,Vsat )

!     write(*,'(I,3F10.6)'),DOY,HOUR,DVS,LAI

     IF(EMR .EQ. .TRUE.)THEN
      !IF(HOUR == 12.5 .AND. DOY < 200)THEN
      !   print *, HGT, LAI, EMR
      !END IF

        CALL CROLST(LLF,LRT,LSTR,WLF,WRT,WSR,WAR,DVS,TRES,hDVS,DLFX1,DLFY1,DLFX2,DLFY2,DLFX3,DLFY3,LAI,LLFst,kLLF,DLFX )

        CALL GRWTH (WAR,GLF,GST,GSR,GSO,GRT,GGR,GWAR,WLF,SWGR,WGR,GRN,GPP-RSP,LSTR,DVS,DVSL,TRES,CFLF,CFST,CFRT,CFSO,hDVS,RTX,RTY,RTX2,LEFY0,LEFX1,LEFY1,LEFX2,LEFY2,LEFX3,LEFY3,PNCLX1,PNCLY1,PNCLX2,PNCLY2,PNCLX3,PNCLY3,FSTR,WSTRS,LAI,LLFst,kLLF,YEAR)

!        write(*,'(2F20.15)'),GLF,WAR

        CALL CROWHT(WLF, WST, WSO, WRT, WSR, WDL,WSH , WAR,WGR , GLF , GST,GSR, GSO, GRT,GGR, LLF, LRT, LSTR,FSTR,CFST)



        CALL CALLAI(LAI,LAImx, WLF, WAR ,DVS,SLWYA,SLWYB,SLWX,CO2PPM)

        CALL CALHGT(HGT,LAI,HGTAA,HGTAB,HGTBA,HGTBB,DVS,hDVS)

        CALL CALROT(ROT,GZRT,MXRT,TRES)

        CALL CALCDI(CDI,TAVE,TMP,DVS,hDVS,PNCLX1,HOUR,TRES,TCmin)

        CALL CALHDI(HDI,TMX,NHED,TMP,DVS,hDVS,HOUR,TRES)

        
     END IF

     

     CALL JUDHVT( YLD,PRM,OCDI,PLT, EMR,GRN, DVS, WLF,WST, WAR ,WSR,WSO,WGR,SWGR, WRT, WDL,NGR,NSP,LAI,HGT,ROT,GDHm1, aGDH,CDI,HDI,NHED,TAVE5,HVT_TAVE,THcrit,HI)

     CALL JUDEMR( WLF, WST,WRT, WAR,DVS, EMR            )

  END IF


  CALL CALSLN(SLN,NFERT,DVS,SLNX1,SLNY1,SLNX2,SLNX3,SLNYMX,SLNYMN,SLNK,CO2PPM)

!  ODVS=DVS


   END SUBROUTINE CROP


!===============================================================
!    Judge plant or not, and initialize variables
!===============================================================


SUBROUTINE JUDPLT (PLT,EMR,GRN,DVS,DVSL,WSH,WLF,WST,WSO,WRT,WSR,WAR,WDL,WGR,SWGR,NGR,NSP,LAI,HGT,ROT,aGDH,aO3FLX,CDI,HDI,NHED,DOY,HOUR,PLANTDOY,IRR,DLFX) !OZN

  IMPLICIT NONE

 !   [OUTPUT]
  INTEGER  PLT         !! Plant or not
  LOGICAL  EMR         !! Emergence or not
  INTEGER  GRN        

  REAL*8   DVS         !! development index
  REAL*8   DVSL

  REAL*8   WSH         !! Weight of shoot

  REAL*8   WLF         !! Weight of leave
  REAL*8   WST         !! Weight of stem
  REAL*8   WSO         !! Weight of storage organ
  REAL*8   WRT         !! Weight of root
  REAL*8   WSR         !! Weight of reserve
  REAL*8   WAR         !! Weight of available reserve
  REAL*8   WDL        !! Weight of dead leave

  REAL*8   WGR
  REAL*8   SWGR

  REAL*8   NGR
  REAL*8   NSP


  REAL*8   LAI
  REAL*8   HGT

  REAL*8   ROT

  REAL*8   aGDH      !! accmulated GDH
  REAL*8   aO3FLX    !! OZN

  REAL*8 CDI
  REAL*8 HDI
  REAL*8 NHED

  REAL*8 DLFX


!   [INPUT]
  INTEGER  DOY         !! # of day from Jan. 1st.
  REAL*8  HOUR
  INTEGER  PLANTDOY         !! planting date as DOY
!  INTEGER  PLTDOY2         !! planting date as DOY

  INTEGER IRR




  IF (    ( DOY .EQ. PLANTDOY )  .and. ( PLT .EQ. 0 )    .and. ( HOUR .GE. 12.D0)  ) THEN

     PLT = 1


     EMR= .FALSE.
     GRN=0
     WLF = 0.D0
     WST = 0.D0
     WRT = 0.D0
     WSO = 0.D0
     WSR = 0.D0
     WAR = 0.D0
     WDL= 0.D0
     WGR = 0.D0
     SWGR = 0.D0

     NGR = 0.D0
     NSP = 0.D0


     WSH = WLF + WST + WSO + WSR + WAR
     aGDH = 0.D0
     aO3FLX = 0.D0
     DVS = 0.D0
     DVS = 0.D0
     ROT = 0.D0

     CDI = 0.D0
     HDI = 0.D0
     NHED = 0.D0

     LAI=0.D0

     DLFX=0.D0
     
  END IF


END SUBROUTINE JUDPLT
!=============================================================
!     Calculation DVS
!=============================================================
!SUBROUTINE CRODVS( DVS   , aGDH, TMP   ,    GDHm1,GDHm2   ,   TRES  ,PLT  )
SUBROUTINE CRODVS( DVS   , aGDH,aVD, TMP   ,    GDHm1 ,   TRES  ,PLT ,VN,TB,TO,TH ,Vsat)
!
  IMPLICIT NONE

!   [MODIFY]
  REAL*8  DVS
  REAL*8  aGDH
  REAL*8 aVD
!
!   [INPUT]
  REAL*8  TMP
  INTEGER  TRES
  REAL*8  GDHm1
!  REAL*8  GDHm2
  INTEGER PLT

  INTEGER VN
  REAL*8 TB
  REAL*8 TO
  REAL*8 TH

  REAL*8 Vsat

!
!   [INTERNAL WORK]
  REAL*8  DVR
!

!   [INTERNAL PARAMETER]
!
  REAL*8 TMPD

  REAL*8,PARAMETER:: TV1=-4.D0
  REAL*8,PARAMETER:: TV2=3.D0
  REAL*8,PARAMETER:: TV3=10.D0
  REAL*8,PARAMETER:: TV4=17.D0

  REAL*8 Veff
  REAL*8 VF
  REAL*8 Vb


  TMPD = TMP - 273.15D0

  IF(TMPD .LT. TB)THEN
     DVR = 0.D0
  ELSE IF(TMPD .LT. TO) THEN
     DVR = TMPD - TB
  ELSE IF(TMPD .LT. TH) THEN
     DVR = (TB-TO)/(TH-TO)*(TMPD-TH)
  ELSE
     DVR = 0.D0
  END IF


  IF(VN>0)THEN
     
     if(TMPD<TV1)THEN
        Veff = 0.D0
     else if(TMPD<TV2)THEN
        Veff = (TMPD-TV1)/(TV2-TV1)
     else if(TMPD<TV3)THEN
        Veff = 1.D0
     else
        Veff = max(0.D0,(TV4-TMPD)/(TV4-TV3))
     END IF
     
     aVD = aVD + Veff * (DBLE(TRES) / 60.D0 / 60.D0) / 24.D0
     
     Vb = Vsat/5.D0
     
     if(aVD<Vb)Then
        VF = 0.D0
     ELSE if(aVD<Vsat)THEN
        VF = (aVD-Vb)/(Vsat-Vb)
     ELSE
        VF = 1.D0
     END IF
     
  ELSE
     VF = 1.0D0
  END IF

  aGDH = aGDH + DVR * (DBLE(TRES) / 60.D0 / 60.D0) / 24.D0 * VF

  DVS = aGDH / GDHm1

END SUBROUTINE CRODVS
!========================================================================
!     Judge  emergence
!========================================================================
SUBROUTINE JUDEMR( WLF, WST, WRT,WAR,DVS, EMR          )
!
  IMPLICIT NONE
!
!    [OUTPUT]
  REAL*8  WLF
  REAL*8  WST
  REAL*8  WRT
  REAL*8  WAR
!
!    [INPUT]
  REAL*8  DVS
  LOGICAL EMR
!
  IF(DVS .GT. 0.012 .AND. EMR .EQ. .FALSE. )THEN
     WLF = 1.D0
     WST = 1.D0
     WRT = 1.D0
     WAR = 0.5D0

     EMR = .TRUE.
  ENDIF

END SUBROUTINE JUDEMR
!========================================================================
!     Loss rate
!========================================================================
SUBROUTINE CROLST( LLF,LRT,LSTR,WLF,WRT,WSR,WAR,DVS,TRES,hDVS,DLFX1,DLFY1,DLFX2,DLFY2,DLFX3,DLFY3,LAI,LLFst,kLLF,DLFX )

  IMPLICIT NONE

  include "./param.inc"
!
!     [OUTPUT]
  REAL*8  LLF
  REAL*8  LRT
  REAL*8  LSTR
!
!     [INPUT]
  REAL*8  WLF
  REAL*8  WRT
  REAL*8  WSR
  REAL*8  WAR
!
  REAL*8  DVS
  INTEGER  TRES

  REAL*8  hDVS
  REAL*8  DLFX1
  REAL*8  DLFY1
  REAL*8  DLFX2
  REAL*8  DLFY2
  REAL*8  DLFX3
  REAL*8  DLFY3

  REAL*8  LAI
  REAL*8  LLFst
  REAL*8  kLLF
  REAL*8  DLFX

!
!    [INTERNAL WORK]
  REAL*8  LLFT
  REAL*8  LSTT
  REAL*8  LRTT
!


!  LLFst=20.D0

!  IF(DLFX==0.D0 .and. LLFst<LAI)THEN
!     DLFX=DVS
!  END IF
!  IF(DLFX==0.D0 .and. hDVS<DVS)THEN
!     DLFX=DVS
!  END IF


!  IF(DLFX>0.D0)THEN
!     LLFT = kLLF * (DVS-DLFX) 
!  ELSE
!     LLFT = 0.D0
!  END IF

  IF(DVS .LT. DLFX1)THEN
     LLFT = DLFY1
  ELSE IF(DVS .LT. DLFX2)THEN
     LLFT = (DLFY2-DLFY1)/(DLFX2-DLFX1)*DVS + DLFY2 - (DLFY2-DLFY1)/(DLFX2-DLFX1)*DLFX2
  ELSE IF(DVS .LT. DLFX3)THEN
     LLFT = (DLFY3-DLFY2)/(DLFX3-DLFX2)*DVS + DLFY3 - (DLFY3-DLFY2)/(DLFX3-DLFX2)*DLFX3
  ELSE
     LLFT = DLFY3
  END IF



!  WRITE(*,'(5F20.10)'),DVS,LAI,LLFT,DLFX,LLFst

  IF(DVS .LT. hDVS)THEN
     LSTT = 0.D0
  ELSE
     LSTT = 10.D0     ! Boumann (2001)      ![day]
  END IF
  LRTT = 0.0

  LLF = (WLF + WAR ) * LLFT * DBLE(TRES)
  LRT = WRT * LRTT * DBLE(TRES)

  IF(LSTT .GT. 0.D0)THEN
     LSTR = WSR * ( 1.D0/ (LSTT * SECPD) ) * DBLE(TRES)
  ELSE
     LSTR=0.D0
  END IF


END SUBROUTINE CROLST


SUBROUTINE GRWTH (  WAR,GLF,GST,GSR,GSO,GRT,GGR,GWAR,WLF,SWGR,WGR,GRN,ASMN,LSTR,DVS,DVSL,TRES,CFLF,CFST,CFRT,CFSO,hDVS,RTX,RTY,RTX2,LEFY0,LEFX1,LEFY1,LEFX2,LEFY2,LEFX3,LEFY3,PNCLX1,PNCLY1,PNCLX2,PNCLY2,PNCLX3,PNCLY3,FSTR,WSTRS ,LAI,LLFst,kLLF,YEAR)
  IMPLICIT NONE

  include "./param.inc"

!   [MODIFY]
  REAL*8  WAR
!
!   [OUTPUT]

  REAL*8  GLF          !! Growth rate of leave            [kg/ha*TRES/SECPD]
  REAL*8  GST          !! Growth rate of stem             [kg/ha*TRES/SECPD]
  REAL*8  GSR
  REAL*8  GSO          !! Growth rate of storage organ    [kg/ha*TRES/SECPD]
  REAL*8  GRT          !! Growth rate of root            [kg/ha*TRES/SECPD]

  REAL*8 GGR      !! Growth rate of grains   [kg/ha*TRES/SECPD]

  REAL*8 GWAR

!
  !   [INPUT]
  REAL*8 WLF

  REAL*8 SWGR
  REAL*8 WGR

  INTEGER GRN

!
  REAL*8 ASMN
  REAL*8 LSTR
!
  REAL*8  DVS
  REAL*8  DVSL
  INTEGER  TRES

  REAL*8  CFLF         ! Fraction for leave
  REAL*8  CFST         ! Fraction for stem
  REAL*8  CFRT         ! Fraction for root
  REAL*8  CFSO         ! Fraction for storage orga

  REAL*8  hDVS
  REAL*8  RTX
  REAL*8  RTY
  REAL*8  RTX2
  REAL*8  LEFY0
  REAL*8  LEFX1
  REAL*8  LEFY1
  REAL*8  LEFX2
  REAL*8  LEFY2
  REAL*8  LEFX3
  REAL*8  LEFY3

  REAL*8  PNCLX1
  REAL*8  PNCLY1
  REAL*8  PNCLX2
  REAL*8  PNCLY2
  REAL*8  PNCLX3
  REAL*8  PNCLY3

  REAL*8 FSTR

  REAL*8  WSTRS

  REAL*8 LAI,kLLF,LLFst

!
!   [INTERNAL WORK]
  REAL*8  CAGCR
  REAL*8  CAGSS
  REAL*8  CAGLF
  REAL*8  CAGST
  REAL*8  CAGSO
  REAL*8  CAGRT
!
  REAL*8  PRTSHT
  REAL*8  PRTLEF
  REAL*8  PRTPNC
!
  REAL*8  GAR          !! Growth rate of reserve            [kg/ha*TRES/SECPD]
!
  REAL*8  REM   !! Remobilization (glucose)

  REAL*8 ppWAR

  REAL*8 LEFXX2
  REAL*8 LEFXX1
  REAL*8 LEFYY1

  INTEGER YEAR

  ppWAR = WAR

  REM = LSTR * CFSR2GL    ! Starch to Glucose           !!! * ST2GLU !* LSSGLU
!
  GAR = ASMN * DBLE(TRES) * 300.D0 ! [mol(CO2)/m**2/s] to [kg(CH2O)/ha/TRES)]

  WAR = WAR + GAR + REM  ! [Glucose]

  WAR = max(WAR,0.D0)

!  write(*,'(2F30.20,F40.30)'),WAR,GAR,ASMN

   IF(WAR .GT. 0.1D0*WLF)THEN  !YM
      CAGCR = (WAR-0.1D0 * WLF) / ((SECPD)*0.5d0/DBLE(TRES))
      !     CAGCR = WAR - 0.1D0 * WLF !YM
   ELSE                         !YM
      CAGCR = 0.D0               !YM
   END IF                       !YM

!!   print *,"CAGCR",CAGCR

  CAGCR = MAX(0.D0,CAGCR)

  WAR = WAR - CAGCR


  GWAR = ppWAR - WAR

  !!! GLUCOSE PARTION !!!
!!  IF(DVS > 0.06D0 .AND. DVS < 0.08D0)THEN
!!     PRTSHT = 0.0D0
!!  ELSE IF(DVS .LT. RTX)THEN
  IF(DVS .LT. RTX)THEN
     PRTSHT = 1.D0 - RTY
     PRTSHT = PRTSHT !* WSTRS * 1.D0
  ELSE IF(DVS .LT. RTX2)THEN
     PRTSHT = 1.D0 - RTY * (RTX2 - DVS)    / (RTX2 - RTX)
     PRTSHT = PRTSHT !* WSTRS * 1.D0
  ELSE
     PRTSHT = 1.D0
  END IF

  




  IF(LAI > LLFst .AND. LEFX1<DVS .AND. GRN==0)THEN
     DVSL = DVS
     GRN = 1
  END IF

  LEFXX2 = (LEFX2 - DVSL)*kLLF + DVSL
  LEFXX1 = DVSL
  LEFYY1 = (LEFY2-LEFY1)/(LEFX2-LEFX1)*DVSL + LEFY2 - (LEFY2-LEFY1)/(LEFX2-LEFX1) * LEFX2

  IF(DVS .LT. LEFX1)THEN
     PRTLEF = (LEFY1-LEFY0)/(LEFX1)*DVS + LEFY0 

  ELSE IF(DVS .LT. LEFX2)THEN

     IF(LAI < LLFst)THEN
        PRTLEF = (LEFY2-LEFY1)/(LEFX2-LEFX1)*DVS + LEFY2 - (LEFY2-LEFY1)/(LEFX2-LEFX1) * LEFX2

     ELSE
        PRTLEF = (LEFY2-LEFYY1)/(LEFXX2-LEFXX1)*DVS + LEFY2 - (LEFY2-LEFYY1)/(LEFXX2-LEFXX1) * LEFXX2
     END IF

  ELSE
     PRTLEF = 0.0D0
  END IF
  PRTLEF=MAX(0.D0,PRTLEF)





!  IF(DVS .LT. LEFX1)THEN
!     PRTLEF = (LEFY1-LEFY0)/(LEFX1)*DVS + LEFY0 

!  ELSE IF(DVS .LT. LEFX2)THEN
!        PRTLEF = (LEFY2-LEFY1)/(LEFX2-LEFX1)*DVS + LEFY2 - (LEFY2-LEFY1)/(LEFX2-LEFX1) * LEFX2
!  ELSE
!     PRTLEF = 0.0D0
!  END IF
!  PRTLEF=MAX(0.D0,PRTLEF)



  IF(DVS .LT. PNCLX1)THEN
    PRTPNC = PNCLY1
  ELSE IF(DVS .LT. PNCLX2)THEN
     PRTPNC = (PNCLY2-PNCLY1)/(PNCLX2-PNCLX1)*DVS + PNCLY2 - (PNCLY2-PNCLY1)/(PNCLX2-PNCLX1) * PNCLX2
  ELSE IF(DVS .LT. PNCLX3)THEN
     PRTPNC = (PNCLY3-PNCLY2)/(PNCLX3-PNCLX2)*DVS + PNCLY3 - (PNCLY3-PNCLY2)/(PNCLX3-PNCLX2) * PNCLX3
  ELSE
     !PRTPNC = 0.5D0 + 0.5D0 / (1960-1896) * (YEAR - 1896)
     !PRTPNC = min(PRTPNC, 1.D0)
     PRTPNC = PNCLY3
  END IF

  IF((PRTLEF+PRTPNC) .GT. 1.D0)THEN    !YM
     PRTLEF = PRTLEF / (PRTLEF+PRTPNC) !YM
     PRTPNC = PRTPNC / (PRTLEF+PRTPNC) !YM
  END IF                               !YM


  CAGSS = CAGCR * PRTSHT                  !!
  CAGRT = CAGCR - CAGSS                  !!
  CAGLF = CAGSS * PRTLEF                   !!
  CAGSO = CAGSS * PRTPNC
  CAGST = CAGSS - CAGLF - CAGSO                 !!
!  CAGST = MIN(CAGST,1.D0)  !YM
!  CAGST = MAX(CAGST,0.D0)  !YM

  GLF = CAGLF * CFLF
!  GST = CAGST * CFST
  GST = CAGST * (1.D0 - FSTR) * CFST
  GSR = CAGST * FSTR / CFSR2GL
  GSO = CAGSO * CFSO
  GRT = CAGRT * CFRT

  IF(DVS > hDVS*0.95D0)THEN
     GGR = GSO
  ELSE
     GGR = 0.D0
  END IF

  !! Sink limitation !!
!  IF(GRN==1)THEN
!  IF(DVS > ((1.D0 - hDVS)*0.2D0 + hDVS))THEN
!     IF(GGR >= (SWGR-WGR))THEN
!        GST = GST + (GGR - (SWGR - WGR))
!        GGR = (SWGR - WGR)
!        GSO = (SWGR - WGR)
!     END IF
!  END IF



END SUBROUTINE GRWTH

!========================================================================
!       Weight crop components
!========================================================================
SUBROUTINE CROWHT(  WLF , WST, WSO, WRT , WSR, WDL ,WSH ,WAR,WGR ,GLF , GST,GSR, GSO,  GRT,GGR ,   LLF , LRT, LSTR,FSTR,CFST  )
!
  IMPLICIT NONE

  include "./param.inc"

!*   [MODIFY]
  REAL*8  WLF       ! Weight of leave                  [kg/ha]
  REAL*8  WST       ! Weight of stem                   [kg/ha]
  REAL*8  WSO       ! Weight of storage organ          [kg/ha]
  REAL*8  WRT       ! Weight of root                   [kg/ha]
  REAL*8  WSR       ! Weight of reserve                [kg/ha]
  REAL*8  WDL      ! Weitht of dead leave
  REAL*8  WGR

!*   [OUTPUT]
  REAL*8  WSH       ! Weight of shoot


!*   [INPUT]
  REAL*8  WAR       !                                  [kg/ha]
  REAL*8  GLF       ! Growth rate of leave             [kg/ha*TRES/SECPD]
  REAL*8  GST       ! Growth rate of stem              [kg/ha*TRES/SECPD]
  REAL*8  GSR
  REAL*8  GSO       ! Growth rate of storage organ     [kg/ha*TRES/SECPD]
  REAL*8  GRT       ! Growth rate of root              [kg/ha*TRES/SECPD]
  REAL*8  GGR

  REAL*8  LLF       ! loss rate of leave               [kg/ha*TRES/SECPD]
  REAL*8  LRT       ! Loss rate of root                [kg/ha*TRES/SECPD]
  REAL*8  LSTR       ! loss rate of sheilded reserve    [kg/ha*TRES/SECPD]



  REAL*8  FSTR
  REAL*8  CFST

  IF( (WLF + GLF - LLF) > 0.D0)THEN
     WLF = WLF + GLF - LLF
  ELSE
     WLF = WLF + GLF
  END IF
!  WST = WST + GST * (1.D0 - FSTR)
  WST = WST + GST 
  WSO = WSO + GSO
  WRT = WRT + GRT - LRT
!  WSR = WSR + GST * FSTR / CFST / CFSR2GL - LSTR
  WSR = WSR + GSR - LSTR

  WGR = WGR + GGR

  WDL = WDL + LLF

  WSH = WLF + WST + WSO  + WSR + WAR

!  print *,"FSTR,WSR",FSTR,WSR,GST

END SUBROUTINE CROWHT
!*============================================================
!*   Calculation LAI at initial and other stages
!*============================================================
SUBROUTINE CALLAI(LAI,LAImx,WLF,WAR,DVS,SLWYA,SLWYB,SLWX,CO2PPM)

  IMPLICIT NONE

  !*   [OUTPUT]
  REAL*8 LAI
  REAL*8 LAImx

  REAL*8 CO2PPM

!*   [INPUT]
  REAL*8 WLF
  REAL*8 WAR
  REAL*8 DVS
  REAL*8 SLWYA
  REAL*8 SLWYB
  REAL*8 SLWX

!*   [INTERNAL WORK]
  REAL*8 SLW

  SLW = SLWYB + (SLWYA - SLWYB) * EXP(-SLWX * DVS) ! org


  SLW = SLW * (1.D0 / (   (0.856D0 * (1.D0 + 1.035D0*exp(-4.35D0*0.001D0*CO2PPM))) / (0.856D0 * (1.D0 + 1.035D0*exp(-4.35D0*0.001D0*368.87D0)))  ) )    !down regulation


  LAI = (WLF + WAR) / SLW

  IF(LAI .GT. 6.D0)THEN
     LAI = 6.D0
  END IF

  IF(LAI .GT. LAImx)THEN
     LAImx = LAI
  END IF



END SUBROUTINE CALLAI
!*====================================================================
!*    Calculation of crop height
!*====================================================================
SUBROUTINE CALHGT( HGT,LAI,HGTAA,HGTAB,HGTBA,HGTBB,DVS,hDVS    )

  IMPLICIT NONE

!*   [OUTPUT]
  REAL*8  HGT

!*   [INPUT]
  REAL*8  LAI
  REAL*8 HGTAA
  REAL*8 HGTAB
  REAL*8 HGTBA
  REAL*8 HGTBB
  REAL*8 DVS
  REAL*8 hDVS

  IF(DVS < hDVS)THEN
     HGT = HGTAA * (DVS / hDVS)
!     HGT = HGTAA * LAI ** HGTBA
  ELSE
     HGT = HGTAA
!     HGT = HGTAB * LAI ** HGTBB
  END IF

!  IF(HGT > 1.5D0)THEN
!     HGT = 1.5D0
!  END IF


END SUBROUTINE CALHGT
!*==============================================================
!*  Calculation root length
!*==============================================================
SUBROUTINE CALROT(ROT,GZRT,MXRT,TRES)

  IMPLICIT NONE

  include "./param.inc"

  REAL*8 GZRT
  REAL*8 MXRT
  INTEGER TRES

  REAL*8 ROT


  ROT = ROT + GZRT / SECPD * DBLE(TRES)
  ROT = MIN(ROT,MXRT)

END SUBROUTINE CALROT
!*===============================================================
!*     Judge harvest or not
!*===============================================================

SUBROUTINE JUDHVT(  YLD,PRM,OCDI,PLT, EMR,GRN, DVS, WLF,WST, WAR,WSR,WSO,WGR,SWGR, WRT,WDL,NGR,NSP,LAI, HGT,ROT,GDHm1, aGDH,CDI,HDI,NHED,TAVE5,HVT_TAVE,THcrit,HI)

  IMPLICIT NONE

!*   [INPUT]
  REAL*8 GDHm1
  REAL*8 aGDH
  REAL*8 CDI
  REAL*8 HDI
  REAL*8 NHED

  REAL*8 TAVE5

  REAL*8 HVT_TAVE

  REAL*8 THcrit

  REAL*8 HI

!*   [OUTPUT]
  REAL*8  YLD
  REAL*8  PRM
  REAL*8  OCDI
  INTEGER PLT
  LOGICAL EMR
  INTEGER GRN
  REAL*8  DVS
  REAL*8  WLF
  REAL*8  WST
  REAL*8  WAR
  REAL*8  WSR
  REAL*8  WSO
  REAL*8  WGR
  REAL*8  SWGR
  REAL*8  WRT
  REAL*8  WDL
  REAL*8  WSH

  REAL*8  NGR
  REAL*8  NSP

  REAL*8  LAI
  REAL*8  HGT
  REAL*8  ROT

  REAL*8 HIC,HIH,HIHC

  

  IF(  ((PLT .EQ. 1)   .AND. (aGDH .GT. GDHm1)) .OR. (  (PLT .EQ. 1) .AND. (TAVE5 < HVT_TAVE)  .AND. (0.5D0<=DVS)   ) ) THEN 


   HIC = HI * (1.D0 - (0.054d0 * (CDI ** 1.56d0))/100.d0)
   IF(HIC < 0.d0)THEN
      HIC = 0.d0
   END IF


   IF(HDI > 0.D0)THEN
      HDI = HDI / NHED - 273.15d0
   ELSE
      HDI = 0.D0
   END IF
      HIH = HI * 1.d0 / (1.d0 + exp(0.853d0 * (HDI - THcrit)))

   IF(HIC < HIH)THEN
      HIHC = HIC
   ELSE
      HIHC = HIH
   END IF

   WSH = WLF + WST + WSO + WSR + WAR
!     YLD1 = WGR / 0.85D0    ! water content 15%
     YLD = WSO * HIHC
     PRM = YLD / WSH
     OCDI = CDI

     DVS = 0.d0
     PLT = 0
     EMR = .FALSE.
     GRN = 0

     WLF = 0.d0
     WST = 0.d0
     WAR = 0.d0
     WSR = 0.d0
     WSO = 0.d0
     WGR = 0.d0
     SWGR = 0.d0
     WRT = 0.d0
     WDL = 0.d0

     NGR = 0.d0
     NSP = 0.d0

     LAI = 0.d0
     HGT = 0.d0
     ROT = 0.D0

     CDI = 0.D0
     HDI = 0.D0
     NHED = 0.D0

     
  END IF
END SUBROUTINE JUDHVT


SUBROUTINE CALCDI(CDI,TAVE,TMP,DVS,hDVS,PNCLX1,HOUR,TRES,TCmin)

  IMPLICIT NONE
  include "./param.inc"

  !*[OUTPUT]
  REAL*8 CDI
  REAL*8 TAVE


  ![INPUT]
  REAL*8 TMP
  REAL*8 DVS
  REAL*8 hDVS
  REAL*8 PNCLX1
  INTEGER TRES
  REAL*8 HOUR

  REAL*8 TCmin


  IF(HOUR < DBLE(TRES)/SECPD*24.d0)THEN
     TAVE = TMP
  ELSE
     TAVE = TAVE + TMP
  END IF


  
  IF( ((PNCLX1+0.05d0) < DVS)  .AND.  (   DVS < ( (1.d0-hDVS)*0.2d0  + hDVS       )  )  )THEN
     IF( (24.d0 - DBLE(TRES)/SECPD*24.d0)    <= HOUR)THEN
        TAVE = TAVE * DBLE(TRES)/SECPD

        IF(TAVE < (TCmin + 273.15d0))THEN
           CDI = CDI + (TCmin + 273.15d0) - TAVE            
        END IF
        
     END IF
  END IF
  
  
END SUBROUTINE CALCDI


SUBROUTINE CALHDI(HDI,TMX,NHED,TMP,DVS,hDVS,HOUR,TRES)



IMPLICIT NONE
include "./param.inc"

!* [OUTPUT]
REAL*8 HDI
REAL*8 TMX
REAL*8 NHED

!* [INPUT]
REAL*8 TMP
REAL*8 DVS
REAL*8 hDVS
REAL*8 HOUR
INTEGER TRES



IF(HOUR < DBLE(TRES)/SECPD*24.d0)THEN
   TMX = TMP
ELSE IF(TMX < TMP)THEN
   TMX = TMP
END IF


IF((hDVS*0.96d0 < DVS) .AND. (   DVS < ( (1.d0-hDVS)*0.2d0  + hDVS       )  )     )THEN

   IF( (24.d0 - DBLE(TRES)/SECPD*24.d0)    <= HOUR)THEN
      HDI = HDI + TMX
      NHED = NHED + 1.d0
   END IF
   
END IF

END SUBROUTINE CALHDI


SUBROUTINE CALSLN(SLN,NFERT,DVS,SLNX1,SLNY1,SLNX2,SLNX3,SLNYMX,SLNYMN,SLNK,CO2PPM)

  IMPLICIT NONE
  include "./param.inc"

!* [OUTPUT]
  REAL*8 SLN

!* [INPUT]
  REAL*8 DVS
  REAL*8 NFERT
  REAL*8 SLNX1
  REAL*8 SLNY1
  REAL*8 SLNX2
  REAL*8 SLNX3
  REAL*8 SLNYMX
  REAL*8 SLNYMN
  REAL*8 SLNK
  
  REAL*8 CO2PPM

  

  REAL*8 X1
  REAL*8 X2
  REAL*8 X3
  REAL*8 Y1
  REAL*8 Y2
  REAL*8 Y3

  


  X1=SLNX1
  X2=SLNX2
  X3=SLNX3
  !Y1=0.7742242627d0
  Y1=SLNY1
!  Y2=3.d0 - (3.d0 - 1.d0) * exp(-0.01d0*NFERT)
  Y2=SLNYMX - (SLNYMX - SLNYMN) * exp(-SLNK*NFERT)
  Y3=0.5D0


  IF(SLNX1 > 0.D0)THEN

     IF(DVS < 0.D0)THEN
        SLN = 0.D0
     ELSE IF(DVS < X1)THEN
        SLN=Y1
     ELSE IF(DVS < X2)THEN
        SLN = (Y1-Y2)*DVS/(X1-X2) + (Y2*X1-Y1*X2)/(X1-X2)
     ELSE IF(DVS < X3)THEN
        SLN = (Y2-Y3)*DVS/(X2-X3) + (Y3*X2-Y2*X3)/(X2-X3)
     ELSE
        SLN = 0.D0
     END IF
     
     
     SLN = SLN * (1.037D0 - 8.33D0 * 0.00001 * CO2PPM)/(1.037D0 - 8.33D0 * 0.00001 * 368.87D0)  !downregulation

  ELSE

     Y3 = (2.25D0 - 1.8D0) / (300.D0) * NFERT + 1.8D0

     IF(DVS < 0.3D0/2.D0)THEN
        SLN = (2.25D0 - 1.5D0)/( 0.3D0/2.D0 ) * (DVS - 0.3D0/2.D0) + 1.5D0
     ELSE IF(DVS < 0.6D0 / 2.D0)THEN
        SLN = (1.7D0 - 2.25D0 ) / (0.6D0/2.D0 - 0.3D0/2.D0) * (DVS - 0.6D0 / 2.D0) + 1.7D0
     ELSE IF(DVS < 1.4D0 / 2.D0)THEN
        SLN = (Y3 - 1.7D0) / (1.4D0/2.D0 - 0.6D0/2.D0) * (DVS - 1.4D0/ 2.D0) + Y3        
     ELSE
        SLN = (0.5d0 - Y3 ) / (1 - 1.4d0 /2.d0) * (DVS - 1.D0) + 0.5D0
     END IF
     
  END IF

!  print *,NFERT,SLNX1,SLNX2,SLNX3,Y2,SLN

END SUBROUTINE CALSLN



