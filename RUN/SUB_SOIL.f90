SUBROUTINE SOIL(WSTRS,WSL,ETC,W2SF,ROT,IRR,GMMSL,TRES,SOILTXT,TMP,PRS,WND,SHM,HGT,PLT,LON,LAT,CRP_NAME) !YM

  IMPLICIT NONE

  include './param.inc'

!  [OUTPUT]
  REAL*8 WSTRS ! water stress for photosynthesis

!  [INPUT]
  REAL*8 ETC ! Transpiration [W/m**2]

  REAL*8 W2SF ! Water to soil surface [kg/m**2/s]
  REAL*8 ROT ! Root length [m]
  INTEGER IRR !
  INTEGER TRES
  REAL*8 BLKS  ! Bulk density of soil
  REAL*8 HCSL  ! Heat conductivity of dry soil
  REAL*8 PORO  ! Effective porosity of soil
  REAL*8 FC    ! field capacity
  REAL*8 WILT  ! Wilting point
  REAL*8 KSAT  !
  REAL*8 PSAT
  REAL*8 BIDX
  REAL*8 GMMSL ! soil water stress factor


  REAL*8 ORGC  ! soil organic content [%]

  INTEGER SOILTXT  !Soil texture
  REAL*8 TMP  ! Tmperature [K]
  REAL*8 PRS ! Air pressure [Pa]
  REAL*8 WND  ! 2m wind [m/s]
  REAL*8 SHM  ! specific humidity
  REAL*8 HGT  ! crop height [m]

  REAL*8 ADSW ! Air-dried soil water: phi=-30000000 [Pa] P107 in Saishindojyogaku
  
  INTEGER PLT  ! 1 or 2 means "plant"

  REAL*8 LON,LAT

  CHARACTER*200 CRP_NAME

  
!  [PARAMETER]
  REAL*8 FG
  REAL*8 WSL0
  DATA FG / 6.D0 /    !MATSIRO defalt
  DATA WSL0 / 0.25D0 /   !MATSIRO defalt

  REAL*8 BRTAU ! time coefficient for Base runoff [day]
  DATA BRTAU / 100.D0 / ! According to Hanasaki et al. (2008)


  REAL*8,PARAMETER::FRACTN=0.02D0  ! Fraction of active humum
  REAL*8,PARAMETER::BETATR=0.00001D0
  REAL*8,PARAMETER::BETAMN=1.4D0
  REAL*8,PARAMETER::BETADENIT=1.4D0
  REAL*8,PARAMETER::GAMSWTH=1.1D0


  !  1: clay(heavy)
  !  2: silty clay
  !  3: clay
  !  4: silty clay loam
  !  5: clay loam
  !  6: silt (=silt loam)
  !  7: silt loam
  !  8: sandy clay
  !  9: loam
  ! 10: sandy clay loam
  ! 11: sandy loam
  ! 12: loamy sand
  ! 13: sand
  REAL*8::PSATS(1:13)=(/-3.7D0,-3.4D0,-3.7D0,-3.3D0,-2.6D0,-2.1D0,-2.1D0,-2.9D0,-1.1D0,-2.8D0,-1.5D0,-0.9D0,-0.7D0/)
  REAL*8::BIDXS(1:13) =(/7.6D0,7.9D0,7.6D0,6.6D0,5.2D0,4.7D0,4.7D0,6.D0,4.5D0,4.D0,3.1D0,2.1D0,1.7D0/)
  REAL*8::KSATS(1:13)=(/0.000017D0,0.000025D0,0.000017D0,0.000042D0,0.000064D0,0.00019D0,0.00019D0,0.000033D0,0.00037D0,0.00012D0,0.00072D0,0.0017D0,0.0058D0/)
  REAL*8::BLKSS(1:13)=(/1330.D0,1260.D0,1330.D0,1300.D0,1390.D0,1380.D0,1380.D0,1470.D0,1430.D0,1500.D0,1460.D0,1430.D0,1430.D0/)
  REAL*8::POROS(1:13) =(/0.5D0,0.52D0,0.5D0,0.51D0,0.48D0,0.48D0,0.48D0,0.44D0,0.46D0,0.43D0,0.45D0,0.46D0,0.46D0/)
  REAL*8::WILTS(1:13) =(/0.3D0,0.27D0,0.3D0,0.22D0,0.22D0,0.06D0,0.11D0,0.25D0,0.14D0,0.17D0,0.08D0,0.05D0,0.05D0/)
  REAL*8::FCS(1:13)   =(/0.42D0,0.41D0,0.42D0,0.38D0,0.36D0,0.3D0,0.31D0,0.36D0,0.28D0,0.27D0,0.18D0,0.12D0,0.1D0/)


!  [INTERNAL VARIABLE]
  INTEGER I,J
  REAL*8 HTCP(NSL)
  REAL*8 HTCD(NSL-1)  ! heat conductivity between soil layers
  REAL*8 HTCDSL(NSL)  ! heat conductivity for all soil layers

  REAL*8 WSL(NSL)
  REAL*8 FAW(NSL)

  REAL*8 Ke !
  REAL*8 HCSLS !

  REAL*8 ETCDFCT !

  REAL*8 FSTRS  !
  REAL*8 EVS ! Evaporation [W/m**2]

  REAL*8 MXEVS
  REAL*8 Rsoil
  REAL*8 Hsoil

  
  REAL*8 WTFLX(0:NSL) ! water flux between each layer; WTFLX(I) indicate water flux from layer I to layer I+1
  REAL*8 K(NSL)   ! water conductivity for each layer
  REAL*8 KB(NSL)  ! water conductivity between layers
  REAL*8 P(NSL)   ! water potential for each layer
  REAL*8 W(NSL)   ! water content ratio to saturation

  REAL*8 DZSL(NSL-1)
  REAL*8 DPTH(0:NSL)

  INTEGER ETLY  ! Index of the lowest layer where root exists
  REAL*8 ETF(NSL)  ! fraction of transpiration for each layer

  ! A*X=B
  REAL*8 A(NSL,NSL)
  REAL*8 B(NSL)
  REAL*8 X(NSL)

  ! Nitrogen dynamics

  REAL*8 ORGCN(NSL) ! Organic C
  REAL*8 NO3(NSL) ! Nitrate concentration
  REAL*8 ORGNH(NSL) ! Organic N (Active + Stable)
  REAL*8 ORGNA(NSL) ! Active Humus
  REAL*8 ORGNS(NSL) ! Stable Humus
  REAL*8 NH4(NSL) ! Anmonium
  SAVE NO3,ORGNA,ORGNS,NH4

  REAL*8 GAMTM(NSL) ! Temperature factor for nutrient cycling
  REAL*8 GAMSW(NSL) ! Soil water factor for nutirent cycling

  REAL*8 NTRANS(NSL) ! The amount of nitrogen transfered from active to stable humus
  REAL*8 NMINE(NSL) ! The amoutn of mineralization
  REAL*8 NDENIT(NSL) ! The amount of denitrification

  REAL*8 NO3PERC(NSL) ! The amont of percolation of NO3



  ! Soil evaporation
  REAL*8 Zd
  REAL*8 Zom
  REAL*8 Zoh
  REAL*8 Ra

  REAL*8 RHO ! air density

  ! [FUNCTION]
  REAL*8 FQSAT


  LOGICAL OFIRST
  SAVE OFIRST
  DATA OFIRST / .TRUE. /




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Thickness between soil layers
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO I=1,NSL-1
     DZSL(I) = (ZSL(I) + ZSL(I+1)) * 0.5D0
  END DO

  DPTH(0)=0.D0
  DO I=1,NSL
     IF(I==1)THEN
        DPTH(I) = ZSL(I)
     ELSE
        DPTH(I) = DPTH(I-1) + ZSL(I)
     END IF
  END DO


  PSAT=PSATS(SOILTXT)
  BIDX=BIDXS(SOILTXT)
  KSAT=KSATS(SOILTXT)
  BLKS=BLKSS(SOILTXT)
  PORO=POROS(SOILTXT)
  WILT=WILTS(SOILTXT)
  FC=FCS(SOILTXT)

  DO I=1,NSL
     WSL(I)=MIN(WSL(I),PORO)
  END DO


  ADSW=PORO*((-30000000.D0/DWATR)/PSAT)**(-1.D0/BIDX)



  !!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!
  !!!!! Soil Water !!!!!
  !!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!


  !!!!!!!!!!!!!!!!!!
  ! Root distribution
  !!!!!!!!!!!!!!!!!!
  IF(ROT>0.D0)THEN
     DO I=1,NSL
        IF(ROT < DPTH(I))THEN
           ETLY = I
           EXIT
        END IF
     END DO


     DO I=1,NSL
        IF(I<ETLY)THEN
           ETF(I)=  3.D0 / 2.D0 / (ROT**3.D0) * (  ((ROT**2.D0) * DPTH(I)-(DPTH(I)**3.D0)/3.D0) - ((ROT**2.D0) * DPTH(I-1)-(DPTH(I-1)**3.D0)/3.D0))
        ELSE IF(I==ETLY)THEN
           ETF(I)=  3.D0 / 2.D0 / (ROT**3.D0) * (  (2.D0 * (ROT**3.D0) /3.D0) - ((ROT**2.D0) * DPTH(I-1)-(DPTH(I-1)**3.D0)/3.D0))
        ELSE
           ETF(I)=0.D0
        END IF
     END DO
  ELSE
     DO I=1,NSL
        ETF(I)=0.D0
     END DO
     ETLY=0
  END IF


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Soil water - Transpiration
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ETCDFCT = 0.D0

!  WRITE(*,'(A10,5F20.10)'),"WSL",((WSL(I)),I=1,5)


  IF(ETLY>0)THEN

  DO I=1,ETLY
     IF( (WSL(I) - WILT) .GT. (ETF(I)*ETC/EL * DBLE(TRES) / DWATR / ZSL(I)))THEN
        WSL(I) = WSL(I)  - ( ETF(I)*ETC/EL * DBLE(TRES) / DWATR / ZSL(I)  )
     ELSE
        ETCDFCT =  ETCDFCT + ETF(I)*ETC  - (WSL(I) - WILT) * EL / DBLE(TRES) * DWATR * ZSL(I)
        WSL(I) = WILT
     END IF
  END DO

!  print *,"ETC",ETC
!  WRITE(*,'(A10,5F20.10)'),"WSL",((WSL(I)),I=1,5)


  IF (ETCDFCT > 0.D0) THEN
     DO I=1,ETLY
        IF( (WSL(I) - WILT) > (ETCDFCT/EL * DBLE(TRES)/DWATR/ZSL(I)) )THEN
           WSL(I) = WSL(I) - (ETCDFCT /EL * DBLE(TRES)/DWATR/ZSL(I))
           ETCDFCT = 0.D0
           EXIT
        ELSE
           ETCDFCT = ETCDFCT - (WSL(I) - WILT)* EL / DBLE(TRES) * DWATR * ZSL(I)
           WSL(I) = WILT
        END IF
     END DO
  END IF

END IF


!  WRITE(*,'(A10,5F20.10)'),"WSL",((WSL(I)),I=1,5)


  EVS=0.D0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!  Evaporation !!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!


  RHO = PRS / RAIR / TMP   ! Air density

  
  
  ! FAO 56 method
  Zd = 2.D0 / 3.D0 * HGT
  Zom = 0.123D0 * HGT
  Zoh = 0.1D0 * Zom



  IF(HGT>0.D0)THEN !YM
     Ra = (log((ZA-Zd)/Zom)*log((ZA-Zd)/Zoh))/(FKARM*FKARM)/WND   ![s/m]
  ELSE             !YM
     Ra = 94.909D0*WND**(-0.9036D0)     !YM  Liu et al. (2007) HESS, Evaluating parameterizations of aerodynamic resistance to heat transfer using field measurements
  END IF           !YM

  !!!!!!!!!!!!!!!!!!!!!!!
  !!! Rsoil and Hsoil !!!
  !!!!!!!!!!!!!!!!!!!!!!!
  
  DO I=1,NSL
     W(I)=WSL(I)/PORO                     ! volumetric water content
     P(I)=PSAT * W(I)**(-BIDX)            ! water potential
  END DO

  
  Rsoil = 800.D0 * (1.D0 - W(1)) / (0.2D0 + W(1))
  Hsoil = exp(P(1)*GRAV/RAIR/TMP)



  !!!! Max EVS  !!!
  MXEVS = WSL(1) * ZSL(1) * DWATR / DBLE(TRES) * EL

!!! Evaporation
  IF(Hsoil*FQSAT(TMP,PRS)>SHM)THEN
     EVS = (Hsoil * FQSAT(TMP,PRS) - SHM)*EL*RHO/(Rsoil+Ra)    ! Sellers 1996   [W/m**2]
     EVS = MIN(EVS,MXEVS)
  ELSE
     EVS = 0.D0
  END IF

!  print *,"EVS",EVS,Ra

  If(IRR==0 .OR. PLT == 0)THEN
     WSL(1) = WSL(1) - EVS / EL * DBLE(TRES) / DWATR / ZSL(1)
  END IF


  DO I=1,NSL
     WSL(I) = MIN(WSL(I),PORO)
     WSL(I) = MAX(WSL(I),ADSW)
  END DO
  
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Calculate coefficients !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO I=1,NSL
     W(I)=WSL(I)/PORO                     ! volumetric water content
     K(I)=KSAT * W(I)**(2.D0*BIDX+3.D0)   ! Hydraulic conductivity
     P(I)=PSAT * W(I)**(-BIDX)            ! water potential

  END DO
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Average Hydraulic conductivity between I and I+1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO I=1,NSL-1
     KB(I)=MIN(K(I),K(I+1))   !YM
  END DO


  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Water flux [kg/m**2/s]
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  BRTAU=100.D0

  DO I=0,NSL
     IF(I==0)THEN
        WTFLX(I)= MIN(W2SF, (PORO-WSL(1))*DWATR*ZSL(1)/DBLE(TRES))
     ELSE IF(I < NSL)THEN
        WTFLX(I)= - KB(I) * ((P(I+1) - P(I)) / DZSL(I) - 1.D0)    ! kumaji
     ELSE
        WTFLX(I)=PORO * DWATR * ZSL(I) / (BRTAU*SECPD) * W(I)**2.D0  !Base runoff based on Hanasaki et al. (2009)
     END IF
  END DO

!  WRITE(*,'(A10,5F20.10)'),"WSL",((WSL(I)),I=1,5)
!  WRITE(*,'(A10,5F20.10)'),"P",((P(I)),I=1,5)
!  WRITE(*,'(A10,6F20.10)'),"WTFLX",((WTFLX(I)),I=0,5)


  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!  Soil water - Flux  !!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO I=1,NSL
     WSL(I) = WSL(I) +  (DBLE(TRES) / ZSL(I) / DWATR) * (WTFLX(I-1) - WTFLX(I))
  END DO



  !!!!!!!!!!!!!!!!
  !! Saturation !!
  !!!!!!!!!!!!!!!!

  IF(IRR==1 .AND. PLT>0)THEN
     DO I=1,NSL
        WSL(I)=PORO
     END DO
  END IF
  DO I=1,NSL
     WSL(I) = MIN(WSL(I),PORO)
     WSL(I) = MAX(WSL(I),ADSW)  !YM

  END DO


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Calculate coefficients !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO I=1,NSL
     W(I)=WSL(I)/PORO                     ! volumetric water content
     P(I)=PSAT * W(I)**(-BIDX)            ! water potential
     FAW(I)=MIN(MAX(WSL(I)-WILT,0.D0)/(FC-WILT),1.D0)  ! Fraction of available water
  END DO


  !!!!!!!!!!!!!!!!!!
  !! Water Stress !!
  !!!!!!!!!!!!!!!!!!
  WSTRS=0.D0
  DO I=1,NSL
     IF( (CRP_NAME .eq. "Rice") .or. (CRP_NAME .eq. "Wheat"))THEN
        FSTRS = 2.D0 / (1.D0 + exp(-P(I)*GMMSL))
     ELSE IF(CRP_NAME .eq. "Soybeans")THEN
        FSTRS = 1.D0 / (1.D0 + 6.24D0*exp(FAW(I)*-12.36D0))
     ELSE
        print *,"Please set Rice, Wheat, or Soybeans"
        stop(0)
     END IF


     WSTRS = WSTRS + FSTRS * ETF(I)
  END DO



END SUBROUTINE SOIL






SUBROUTINE GAUSS_JORDAN(A0,B,X,N)

  IMPLICIT NONE

!  [OUTPUT]
  REAL*8 X(N)

!  [INPUT]
  INTEGER N
  REAL*8 A0(N,N)
  REAL*8 B(N)

!  [INTERNAL VARIABLE]
  INTEGER I,K
  REAL*8 AR
  REAL*8 A(N,N)

  A(:,:)=A0(:,:)
  X(:)=B(:)

  DO K = 1,N
     IF(A(K,K)==0.0D0) STOP
     AR = 1.0D0 / A(K,K)
     A(K,K) = 1.0D0
     A(K,K+1:N)=AR*A(K,K+1:N)
     X(K) = AR * X(K)
     DO I = 1,N
        IF(I /= K) THEN
           A(I,K+1:N) = A(I,K+1:N) - A(I,K) * A(K,K+1:N)
           X(I) = X(I) - A(I,K) * X(K)
           A(I,K) = 0.0D0
        END IF
     END DO
  END DO

END SUBROUTINE GAUSS_JORDAN
