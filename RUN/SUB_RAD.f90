SUBROUTINE RAD(QPARSNLF,QPARSHLF,VMXSNLF,VMXSHLF,LAISN,LAISH,SLN,KN,RSD,LAI,RLFV,TLFV,RLFN,TLFN,LAT,DOY,HOUR,ILON,ILAT,CRP_NAME) !YM
  
  IMPLICIT NONE

  include "./param.inc"

!  [OUTPUT]
  REAL*8 QPARSNLF   ! Absorbed PAR by sunlit leaf per unit leaf area [W/m**2(l)]
  REAL*8 QPARSHLF   ! Absorbed PAR by shade leaf per unit leaf area [W/m**2(l)]
  REAL*8 VMXSNLF
  REAL*8 VMXSHLF
  REAL*8 LAISN
  REAL*8 LAISH


!  [INPUT]
  REAL*8 SLN
  REAL*8 VMXtop
  REAL*8 KN
  REAL*8 PAR  ! PAR                                             [W/m**2]
  REAL*8 RSD  ! Downward short wave radiation                   [W/m**2]
  REAL*8 LAI  ! Total leaf area (Sunlit + Shaded)               [m**2(l)/m**2]
  REAL*8 RLFV,TLFV ! Refrection and Transmittion coefficients of leaf for PAR [-]
  REAL*8 RLFN,TLFN ! Refrection and Transmittion coefficients of leaf for Near Infrared [-]
  REAL*8 LAT       ! Latitude [degree] (**.** (-90.0 - 90.0))
  INTEGER DOY      ! Day of year [Day]
  REAL*8 HOUR      ! Hour        [Hour]

  INTEGER ILAT,ILON
  
  CHARACTER*200 CRP_NAME

!  [INTERNAL VARIABLE]
  INTEGER I
  REAL*8 COSP
  REAL*8 SECP

  REAL*8 RLF(2)     !leaf albedo 1:PAR,2:NIR
  REAL*8 TLF(2)    !leaf trans. 1:PAR,2:NIR


  REAL*8 NIR
  REAL*8 DIR(2)   ! 1:PAR,2:NIR
  REAL*8 DIF(2)   ! 1:PAR,2:NIR

  REAL*8 SC
  REAL*8 ATMTR
  REAL*8 FRDF

  REAL*8 EAP(2)
  REAL*8 EAN(2)
  REAL*8 SRF(2)
  REAL*8 RFE(2)
  REAL*8 TFE(2)

  REAL*8 DIRLAI
  REAL*8 a(2)
  REAL*8 A1(2)
  REAL*8 A2(2)
  REAL*8 A3(2)
  REAL*8 C1(2)
  REAL*8 C2(2)
  REAL*8 C3(2)
  REAL*8 C4(2)
  
  REAL*8 ALBCDIF(2)
  REAL*8 ALBCDIR(2)
  REAL*8 TRCDIF(2)
  REAL*8 TRCDIR(2)

  REAL*8 QSN
  REAL*8 QSH
  REAL*8 QSNDIR
  REAL*8 QSNDIF
  REAL*8 QSHDIF
  REAL*8 QDIF

  REAL*8 VMX
  REAL*8 VMXSN
  REAL*8 VMXSH

!  [FUNCTION]
  REAL*8 SINB

!  [INTERNAL PARAM]
  REAL*8 DF     ! Diffusivity factor
  REAL*8 FLFE   ! effective leaf param



  DATA DF / 1.66D0 /    ! Acording to Kondoh (1994) P70
  DATA FLFE / 1.D0 /    ! FLFE is 0.75 in the Original MATSIRO 


  !!!!!!!!!!!!!!!!!!!!
  !!! TRCL 
  !!!!!!!!!!!!!!!!!!!!
!  TRCL = EXP(-LFOR * LAI * DF)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Solar elevation (beta) 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  COSP = SINB(DOY,HOUR,LAT)    ! Cos(phi) = Sin(beta) , phi = pi/2 - beta
!  IF(COSP==0.D0)THEN
!     print *,"cosp"
!     stop
!  END IF
!YM  SECP = 1.D0 / COSP


!  IF( COSP .LE. 0.D0) THEN
!     RSD = 0.D0
!     PAR = 0.D0
!  END IF


  TLF(1) = TLFV
  RLF(1) = RLFV
  TLF(2) = TLFN
  RLF(2) = RLFN


  IF( LAI .GT. 0.D0)THEN

     DO I=1,2  ! PAR:1,NIR:2
        RFE(I) = FLFE * RLF(I) + (1.D0 - FLFE) * TLF(I)
        TFE(I) = FLFE * TLF(I) + (1.D0 - FLFE) * RLF(I)

        SRF(I) = SQRT( (1.D0-TFE(I))**2.D0 - RFE(I)**2.D0)
        a(I) = LFOR*DF*SRF(I)
  
     END DO

     IF( (CRP_NAME .eq. "Rice") .or. (CRP_NAME .eq. "Wheat"))THEN

        VMXtop = MIN(87.04D0*(SLN-0.487D0),138.77D0)  ! Rice

     ELSE IF(CRP_NAME .eq. "Soybeans")THEN

        VMXtop = MIN(70.D0*(SLN-0.5D0),140.D0)  ! Soy

     ELSE
        
        print *,"Please set Rice, Wheat, or Soybeans"
        stop(0)

     END IF
!     VMXtop = MIN(87.04D0*(2.D0/3.D0)*(SLN-0.487D0),138.77D0*(2.D0/3.D0))  ! Soybean




!     VMXtop = 100.D0
!     print *,"VMXtop",VMXtop

        VMX = (VMXtop/1000000.D0) * (1.D0 - EXP(-KN*LAI))/KN   ! CLM
!        VMX = (VMXtop/1000000.D0) * (1.D0 - EXP(-KN))/KN * LAI  ! de Pury et al.

!        IF((ILAT==20).and.(ILON==33))THEN
!           print *,VMX,VMXtop,LAI,KN
!        END IF

           
!        IF( COSP .GT. 0.D0 )THEN
!        IF( COSP .GT. 0.D0 .AND. ( (RSD .GT. 0.D0)  .OR. (PAR .GT. 0.D0) ))THEN    ! modifiedin 26 May, 2018
         IF(   ( (RSD .GT. 0.D0)  .OR. (PAR .GT. 0.D0) ))THEN    ! modifiedin 4 Sep, 2019
        
!        SECP = 1.D0 / COSP !YM

        
        !!!!!!!!!!!!!
        ! NIR
        !!!!!!!!!!!!!
!
        PAR = 0.5D0 * RSD
        NIR = RSD - PAR
        IF(NIR .LT. 0.D0)THEN
           NIR = 0.D0
           RSD = PAR
        END IF


        !*******************************************************
        ! Direct and Diffuse by Goudriaan and van Laar (1994)
        !*******************************************************
        SC = 1370.D0 * (1.D0 + 0.033D0 * COS(2.0*PI*DBLE(DOY)/365.D0))  !Solar constant


        IF(COSP .GT. 0.D0)THEN
           SECP = 1.D0 / COSP
           ATMTR = RSD / (SC*COSP)
           IF(ATMTR .LE. 0.22D0)THEN
              FRDF = 1.D0
           ELSE IF (ATMTR .LE. 0.35D0)THEN
              FRDF = 1.D0 - 6.4D0 * (ATMTR - 0.22D0) ** 2.D0
           ELSE
              FRDF = 1.47D0 - 1.66D0 * ATMTR
           END IF
           FRDF = MAX(FRDF, 0.15D0 + 0.85D0 * (1.D0 - EXP(-0.1D0 / COSP)))


        !!!!!!!!!!!!!!!!!!!!!!
        !!!!! Global MATCRO for S12
!        IF(RSD > 0.D0)THEN
!           FRDF = RSDDIF/RSD
!        ELSE
!           FRDF = 0.D0
!        END IF
        !!!!!!!!!!

           DIF(1)= PAR * FRDF
           DIR(1)= PAR - DIF(1)
           DIF(2)= NIR * FRDF
           DIR(2)= NIR - DIF(2)



           DO I=1,2  ! PAR:1,NIR:2



              EAP(I) = EXP(a(I)*LAI)
              EAN(I) = EXP(-a(I)*LAI)

              DIRLAI = DIR(I) * EXP(-LFOR * LAI * SECP)

              A1(I) = ( 1.D0 - TFE(I) + SRF(I)  )/RFE(I)
              A2(I) = ( 1.D0 - TFE(I) - SRF(I)  )/RFE(I)
              A3(I) = (A1(I) - ALBS) * EAP(I) - (A2(I) - ALBS) * EAN(I)
              C3(I) = SECP * (TFE(I)*SECP+DF*TFE(I)*(1.D0-TFE(I)) + DF * RFE(I)**2.D0   ) / (DF**2.D0 * ( (1.D0 - TFE(I))**2.D0 - RFE(I)**2.D0  ) -SECP**2.D0 )
              C4(I) = (RFE(I) * (DF - SECP) * SECP) / (DF**2.D0 * ( (1.D0 - TFE(I))**2.D0 - RFE(I)**2.D0  ) -SECP**2.D0 )

              C1(I) = ( -(A2(I) - ALBS)*EAN(I)*(DIF(I)-C3(I)*DIR(I)) + (C3(I)*ALBS+ALBS-C4(I))*DIRLAI)/A3(I)
              C2(I) = (  (A1(I) - ALBS)*EAP(I)*(DIF(I)-C3(I)*DIR(I)) - (C3(I)*ALBS+ALBS-C4(I))*DIRLAI)/A3(I)
              
           ENDDO
  

        !***************************************
        ! LAI in Sunlit(LAISN) and Shade(LAISH)
        !***************************************
           LAISN = (1.D0 -EXP(-LFOR*SECP*LAI)) / (LFOR*SECP)
           LAISH = LAI - LAISN




        !**************************************************
        ! Absorbed PAR  per unit leaf area 
        !**************************************************

           QSNDIR = DIR(1) * (1.D0 - EXP(-LFOR*SECP*LAI))  

           QDIF = C1(1)*(1.D0 -A1(1))*(1.D0-EXP(a(1)*LAI)) + C2(1)*(1.D0-A2(1))*(1.D0-EXP(-a(1)*LAI)) +(C3(1)-C4(1))*DIR(1)*(1.D0-EXP(-LFOR*SECP*LAI))
           QDIF = MAX(QDIF,0.D0)

           QSNDIF = a(1)*C1(1)*(1.D0-A1(1))/(a(1)-LFOR*SECP)*(1.D0-EXP((a(1)-LFOR*SECP)*LAI)) + a(1)*C2(1)*(1.D0-A2(1))/(a(1)+LFOR*SECP)*(1.D0-EXP(-(a(1)+LFOR*SECP)*LAI)) + (C3(1)-C4(1))*DIR(1)/2.D0*(1.D0-EXP(-2.D0*LFOR*SECP*LAI))
           QSNDIF = MAX(QSNDIF,0.D0)

!        print *,"QDIF",QDIF,QSNDIF

           QSHDIF = QDIF - QSNDIF
           QSHDIF = MAX(QSHDIF,0.D0)


           QSN = QSNDIR + QSNDIF
           QSH = QSHDIF


           IF(LAISN > 0.D0)THEN
              QPARSNLF = QSN / LAISN
           ELSE
              QPARSNLF = 0.D0
           END IF

           IF(LAISH > 0.D0)THEN
              QPARSHLF = QSH / LAISH
           ELSE
              QPARSHLF = 0.D0
           END IF






        !****************************************
        ! Vmax in Sunlit(VMXSN) and Shade(VMNSH)
        !****************************************

           VMXSN = (VMXtop/1000000.D0) * (1.D0 - EXP(-(KN+LFOR*SECP)*LAI)) / (KN+LFOR*SECP)  ! CLM
!        VMXSN = (VMXtop/1000000.D0) * (1.D0 - EXP(-(KN+LFOR*SECP*LAI))) / (KN+LFOR*SECP*LAI) * LAI  ! de Pury 
           VMXSN = MIN(VMX,VMXSN)


           VMXSH = VMX - VMXSN


           IF(LAISN > 0.D0)THEN
              VMXSNLF = VMXSN / LAISN
           ELSE
              VMXSNLF = 0.D0
           END IF

           IF(LAISH > 0.D0)THEN
              VMXSHLF = VMXSH / LAISH
           ELSE
              VMXSHLF = 0.D0
           END IF

        ELSE

           FRDF = 1.D0

           DIF(1)= PAR * FRDF
           DIR(1)= PAR - DIF(1)
           DIF(2)= NIR * FRDF
           DIR(2)= NIR - DIF(2)



           DO I=1,2  ! PAR:1,NIR:2


              EAP(I) = EXP(a(I)*LAI)
              EAN(I) = EXP(-a(I)*LAI)


              A1(I) = ( 1.D0 - TFE(I) + SRF(I)  )/RFE(I)
              A2(I) = ( 1.D0 - TFE(I) - SRF(I)  )/RFE(I)
              A3(I) = (A1(I) - ALBS) * EAP(I) - (A2(I) - ALBS) * EAN(I)


              C1(I) = ( -(A2(I) - ALBS)*EAN(I)*DIF(I))/A3(I)
              C2(I) = (  (A1(I) - ALBS)*EAP(I)*DIF(I))/A3(I)
              
           ENDDO
  


        !**************************************************
        ! Absorbed PAR  per unit leaf area 
        !**************************************************


           QDIF = C1(1)*(1.D0 -A1(1))*(1.D0-EXP(a(1)*LAI)) + C2(1)*(1.D0-A2(1))*(1.D0-EXP(-a(1)*LAI))
           QDIF = MAX(QDIF,0.D0)

           QSNDIF = 0.D0

           QSHDIF = QDIF - QSNDIF
           QSHDIF = MAX(QSHDIF,0.D0)
           QSH = QSHDIF


           LAISN = 0.D0
           LAISH = LAI

           QPARSNLF = 0.D0
           QPARSHLF = QSH / LAISH
           VMXSNLF = 0.D0
           VMXSHLF = VMX / LAISH

        END IF


     ELSE

        QPARSNLF=0.D0
        QPARSHLF=0.D0
        LAISN = 0.D0
        LAISH = LAI
        VMXSNLF = 0.D0
        VMXSHLF = VMX / LAI

     END IF


  ELSE

     QPARSNLF=0.D0
     QPARSHLF=0.D0
     VMXSNLF = 0.D0
     VMXSHLF = 0.D0
     LAISN = 0.D0
     LAISH = 0.D0

  END IF


END SUBROUTINE RAD
