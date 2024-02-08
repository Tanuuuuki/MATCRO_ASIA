PROGRAM MATCRO

  IMPLICIT NONE
  include "./param.inc"

  !! !!!!!!!!!!!!!!!!!!!!!!! !!
  !! [[ SETTING PARAMETER ]] !!
  !! !!!!!!!!!!!!!!!!!!!!!!! !!
  

  !! [TIME and REGION] !!
  CHARACTER*5 MODE  ! MODE for Calculation : POINT or SURFACE
  CHARACTER*5 POINT ! Calculating point [-]
  INTEGER STYR      ! Start year [Year]
  INTEGER ENYR      ! End year   [Year]
  INTEGER STDOY     ! Start DOY  [day]
  INTEGER ENDOY     ! End DOY    [day]
  

  INTEGER TRES      ! Time resolution     [Second]
  INTEGER clmTRES  ! Time resolution of climate forcing [Second]
  


  REAL*8  WEST      ! West edge  [degree]
  REAL*8  EAST      ! East edge  [degree]
  REAL*8  NORTH     ! North edge [degree]
  REAL*8  SOUTH     ! South edge [degree]
  REAL*8  WERES     ! Resolution of west to east [degree]
  REAL*8  NSRES     ! Resolution of north to south [degree]
  

  !! [CROP NAME and FILE] !!
  CHARACTER*200 CRP_NAME   ! Crop name (Rice, Wheat, Maize, or Soybeans)
  CHARACTER*200 CRP_FILE   ! Input file name for crop parameter
  
  !! [CLIMATE FILE] !!
  CHARACTER*200 CLIM_HEAD  ! Head of input file name for climate data
  CHARACTER*200 CLIM_FOOT  ! Foot of input file name for climate data  
  
  !! [INPUT FILE] !!
  CHARACTER*200 CNFILE     ! Head of input file name for atmospheric CO2
  CHARACTER*200 PLT_FILE   ! Foot of input file name for atmospheric CO2
  CHARACTER*200 SMPL_FILE  ! Foot of input file name for atmospheric CO2
  
  !! [OUTPUT FILE] !!
  CHARACTER*200 YLD_FILE   ! Output file name for yields of first growing season
  CHARACTER*200 PRM_FILE   ! Output file name for PRM ← HI for final
  CHARACTER*200 CDI_FILE   ! Output file name for PRM ← HI for final
  
  !! [OTHERs] !!
  REAL*8  WND_HGT          ! Hegith of input data for input [m]
  !INTEGER NGRW             ! Number of growing season in a year
  INTEGER IER
  
  !! !!!!!!!!!!!!!!!!!!! !!
  !! [LOCATION VARIABLE] !!
  !! !!!!!!!!!!!!!!!!!!! !!
  REAL*8 LON  ! Longitude
  REAL*8 LAT  ! Latitude
  REAL*8 POINTLON   ! Longitude for point calculation
  REAL*8 POINTLAT   ! Latitude for point calculation
  INTEGER ILAT  ! Index of latitude
  INTEGER ILON  ! Index of longitude
  INTEGER NLAT  ! Number of grids in latitude
  INTEGER NLON  ! Number of grids in longitude
  

  !! !!!!!!!!!!!!!!!!!! !!
  !! [FORCING VARIABLE] !!
  !! !!!!!!!!!!!!!!!!!! !!
  REAL*8 TMP    ! Air temperature              [K]
  REAL*8 TMX    ! Daily Maximum temperature    [K]
  REAL*8 PRC    ! Precipitation                [kg/m**2/s]
  REAL*8 RSD    ! Downward shortwave radiation [W/m-2]
  REAL*8 SHM    ! Specific density of H2O      [kg/kg]
  REAL*8 WND    ! Wind speed                   [m/s]
  REAL*8 PRS    ! Surface pressure             [Pa]
  REAL*8 OZN    ! Ozone conentration           [ppb]

  REAL*8,ALLOCATABLE:: INTMP(:,:,:)    ! Air temperature               [same unit of the input data]
  REAL*8,ALLOCATABLE:: INTMX(:,:,:)    ! Maximum Air temperature       [same unit of the input data]
  REAL*8,ALLOCATABLE:: INTMN(:,:,:)    ! Minimum Air temperature       [same unit of the input data]
  REAL*8,ALLOCATABLE:: INPRC(:,:,:)    ! Precipitation                 [same unit of the input data]
  REAL*8,ALLOCATABLE:: INRSD(:,:,:)    ! Downward shortwave radiation  [same unit of the input data]
  REAL*8,ALLOCATABLE:: INSHM(:,:,:)    ! Specific density of H2O       [same unit of the input data]
  REAL*8,ALLOCATABLE:: INWND(:,:,:)    ! Wind speed                    [same unit of the input data]
  REAL*8,ALLOCATABLE:: INPRS(:,:,:)    ! Surface pressure              [same unit of the input data]
  REAL*8,ALLOCATABLE:: INOZN(:,:,:)    ! Ozone concentration           [same unit of the input data]



  REAL*8 CO2PPM ! CO2 concentration [PPM]

  !! !!!!!!!!!!!!!!!! !!
  !! [CROP PARAMETER] !!
  !! !!!!!!!!!!!!!!!! !!
  !<Leaf photosynthesis>
  REAL*8 RESPCP    ! Respiration fraction of Vmax                   [-]
  REAL*8 EFFCON    ! Quantum efficiency                             [mol/mol]
  REAL*8 ATHETA    ! Coupling parameter (Wc,We)                     [-]
  REAL*8 BTHETA    ! Coupling parameter (Wp,Ws)                     [-]
  REAL*8 MH2O      ! Conductance-phtosynthesis slope parameter      [mol/m**2(l)/s] (for GSH2O)
  REAL*8 BH2O      ! Conductance-photosynthesis intercept           [mol/m**2(l)/s] (for GSH2O)
  REAL*8 KN        ! Vertical distribution factor for Nitrogen      [-]

  REAL*8 ZKCA      ! Kc at 298K                                     [Pa]
  REAL*8 ZKCB      ! Parameter for temperature dependence of Kc     [-]
  REAL*8 ZKOA      ! Ko at 298K                                     [Pa]
  REAL*8 ZKOB      ! Parameter nfor temperature dependence of Ko    [-]
  REAL*8 GMMA      ! Parameter for temperature dependence of Gamma* [-]
  REAL*8 GMMB      ! Parameter for temperature dependence of Gamma* [-]

  !<Bulk transfer coefficient>
  REAL*8 LTCH      ! Leaf transfer coefficient for heat             [-]

  !<Radiation>
  REAL*8 RLFV      ! Leaf albedo (VIS)                              [-]
  REAL*8 TLFV      ! Leaf trans. (VIS)                              [-]
  REAL*8 RLFN      ! Leaf albedo (NIR)                              [-]
  REAL*8 TLFN      ! Leaf trans. (NIR)                              [-]

  !<Crop>
  REAL*8 hDVS      ! DVS at heading                                 [-]

  REAL*8 CFLF      ! Fraction of C for leaf                         [kg(C)/kg(leaf)]
  REAL*8 CFST      ! Fraction of C for stem                         [kg(C)/kg(stem)]
  REAL*8 CFRT      ! Fraction of C for root                         [kg(C)/kg(root)]
  REAL*8 CFSO      ! Fraction of C for storage                      [kg(C)/kg(storage)]

  INTEGER VN
  REAL*8 TB
  REAL*8 TO
  REAL*8 TH

  REAL*8 DLFX1     ! 1st point of DVS for dead leaf                 [-]
  REAL*8 DLFY1     ! Rate of dead leaf at 1st point of DVS          [-]
  REAL*8 DLFX2     ! 2nd point of DVS for dead leaf                 [-]
  REAL*8 DLFY2     ! Rate of dead leaf at 2nd point of DVS          [-]
  REAL*8 DLFX3     ! 3rd point of DVS for dead leaf                 [-]
  REAL*8 DLFY3     ! Rate of dead leaf at 3rd point of DVS          [-]

  REAL*8 RTX       ! 1st point of DVS for partition for root        [-]
  REAL*8 RTY       ! Rate of partition for root at 1st point        [-]
  REAL*8 RTX2      ! 2nd point of DVS for partition for root        [-]

  REAL*8 LEFY0     ! Rate of partition for leaf at DVS=0            [-]
  REAL*8 LEFX1     ! 1st point of DVS for leaf                      [-]
  REAL*8 LEFY1     ! Rate of partition for leaf at 1st point        [-]
  REAL*8 LEFX2     ! 2nd point of DVS for leaf                      [-]
  REAL*8 LEFY2     ! Rate of partition for leaf at 2nd point        [-]
  REAL*8 LEFX3     ! 3rd point of DVS for leaf                      [-]
  REAL*8 LEFY3     ! Rate of partition ofr leaf at 3rd point        [-]

  REAL*8 PNCLX1    ! 1st point of DVS for panicle                   [-]
  REAL*8 PNCLY1    ! Rate of partition for panicle at 1st point     [-]
  REAL*8 PNCLX2    ! 2nd point of DVS for panicle                   [-]
  REAL*8 PNCLY2    ! Rate of partition for panicle at 2nd point     [-]
  REAL*8 PNCLX3    ! 3rd point of DVS for panicle                   [-]
  REAL*8 PNCLY3    ! Rate of partition for panicle at 3rd point     [-]

  REAL*8 FSTR      ! Partition to sielded reserve (starch) in stem  [-]

  REAL*8 SLWYA     ! Specific leaf weight at DVS = 0                [kg/m**2(l)]
  REAL*8 SLWYB     ! Specific leaf weight at  DVS = infinity        [kg/m**2(l)]
  REAL*8 SLWX      ! Shape parameter of specific leaf weight        [-]

  REAL*8 HGTAA     ! HGTAA * LAI ** HGTBA (before heading)          [-]
  REAL*8 HGTAB     ! HGTAB * LAI ** HGTBB (after heading)           [-]
  REAL*8 HGTBA     ! HGTAA * LAI ** HGTBA (before heading)          [-]
  REAL*8 HGTBB     ! HGTAB * LAI ** HGTBB (after heading)           [-]

  REAL*8 GZRT      ! Growth rate of root                            [m/day]
  REAL*8 MXRT      ! Maximum root length                            [m]
  REAL*8 GMMSL     ! Soil stress factor                             [-]

  REAL*8 SLNX1     ! 1st point of DVS for SLN                       [-]
  REAL*8 SLNY1     ! 1st point of DVS for SLN                       [-]
  REAL*8 SLNX2     ! 2nd point of DVS for SLN                       [-]
  REAL*8 SLNX3     ! 3rd point of DVS for SLN                       [-]
  REAL*8 SLNYMX    ! Maximum SLN                                    [-]
  REAL*8 SLNYMN    ! Minimum SLN                                    [-]
  REAL*8 SLNK      ! Shape parameter of SLN curve                   [-]

  REAL*8 TCmin     ! Critical temperature for cool damage           [deg] 
  REAL*8 THcrit    ! Critical temperature for hot damage            [deg]
  REAL*8 HI
 
  REAL*8 PLTDIF    ! Day of difference of PLTDOY                    [Day]

  REAL*8 HVT_TAVE  ! Critical temperature at harvest 

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!
  !! [PLANT DOY and GDH at maturity] !!
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!
  INTEGER,ALLOCATABLE:: SOILTXT(:,:)    ! Soil texture            [-]
  INTEGER,ALLOCATABLE:: PLTDOY(:,:)   ! Planting doy            [DOY]
  REAL*8, ALLOCATABLE:: GDHm(:,:)     ! GDH at maturity         [Degree seconds]
  REAL*8, ALLOCATABLE:: GDHh(:,:)     ! GDH at heading          [Degree seconds]
  REAL*8, ALLOCATABLE:: NFERT(:,:)      ! Nitrogen fertilizer application [kg/ha]
  REAL*8, ALLOCATABLE:: Vsat(:,:)        ! required vernarization days
  REAL*8, ALLOCATABLE:: LND(:,:)        ! required vernarization days

  REAL*8, ALLOCATABLE:: BUF(:,:)        ! Temerate array for inputs [same unit of inputs]

  !! !!!!!!!!!!!!!!!!!!! !!
  !! [INTERNAL VARIABLE] !!
  !! !!!!!!!!!!!!!!!!!!! !!
  INTEGER IRR         ! FLAG of irrigation     [-] IRR=1: IRRIGATION; IRR=0: RAINFED
  !INTEGER GRW         ! FLAG of growing season [-] GRW=1: first growing season; GRW=2: second growing season  
  !INTEGER IGRW
  

  INTEGER YEAR        ! Year                   [Year]
  INTEGER DOY         ! Day of year            [Day]
  REAL*8  HOUR        ! Hour                   [Hour]

  INTEGER IHOUR       ! Index of hour          [-]
  INTEGER NDOY        ! DOY of next day        [-] 
  INTEGER PDOY        ! DOY of previous day    [-]

  REAL*8 INTSINB      ! Integration of SINB during a day [second]


  REAL*8 GPP          ! Gross primary production         [mol(CO2)/m**2/s]
  REAL*8 RSP          ! Respiration                      [mol(CO2)/m**2/s]
  REAL*8 TSP          ! Transpiration                    [kg/m**2/s]

  REAL*8 QPARSNLF     ! Absorbed PAR per unit leaf area for sunlit leaf   [W/m**2(l)]
  REAL*8 QPARSHLF     ! Absorbed PAR per unit leaf area for sunshade leaf [W/m**2(l)]
  REAL*8 VMXSNLF      ! Vmax per unit leaf area for sunlit leaf           [CO2/m**2(l)/s]
  REAL*8 VMXSHLF      ! Vmax per unit leaf area for sunshade leaf         [CO2/m**2(l)/s]
  REAL*8 LAISN        ! LAI for sunlit leaf                               [m**2(l)/m**2]
  REAL*8 LAISH        ! LAI for sunshade leaf                             [m**2(l)/m**2]

  INTEGER I,J

  INTEGER STDD,ENDD

  !! < State variable > !!
  REAL*8 WSL(NSL)     ! Soil water content       [m**3/m**3] 
  REAL*8 WSTRS        ! Soil water stress factor [-] (0~1)

  REAL*8 SLN          ! Leaf nitrogen concentration per unit leaf area[g/m**2(l)]

  INTEGER PLT         ! FLAG for planting        [-]
  LOGICAL EMR         ! FLAG for emergence       [-]
  INTEGER GRN         ! FLAG for grain filling   [-]

  REAL*8 DVS          ! OZN
  REAL*8 DVSL
  REAL*8 aO3FLX       ! OZN
  REAL*8 aGDH         ! accumulated GDH

  REAL*8 aVD  

  REAL*8 CDI          ! Cold damage index
  REAL*8 HDI          ! Hot damage index
  REAL*8 TAVE
  REAL*8 NHED         ! Number of days during heading for heat damage [day]


  REAL*8 LAI          ! LAI               [m**2(l)/m**2(ground)]
  REAL*8 ROT          ! Root depth             [m]
  REAL*8 HGT          ! Height            [m]

  REAL*8 DLFX
  REAL*8 LLFst
  REAL*8 kLLF

  REAL*8 WSH          ! Weight of shoot                                 [kg/ha]
  REAL*8 WSO          ! Weight of storage organ                         [kg/ha]
  REAL*8 WST          ! Weight of stem                                  [kg/ha]
  REAL*8 WLF          ! Weight of leaf                                  [kg/ha]
  REAL*8 WRT          ! Weight of root                                  [kg/ha]
  REAL*8 WAR          ! Weight of easily available reserve (glucose)    [kg/ha]
  REAL*8 WIR          ! Weight of non-easily available reserve (starch) [kg/ha]
  REAL*8 WDL          ! Weight of dead leaf                             [kg/ha]
  REAL*8 WGR
  REAL*8 SWGR
  REAL*8 NGR
  REAL*8 NSP


  REAL*8 GLF
  REAL*8 GST
  REAL*8 GRT
  REAL*8 GSO
  

! Storage for state variables
  REAL*8,ALLOCATABLE:: pWSL(:,:,:)
  REAL*8,ALLOCATABLE:: pWSTRS(:,:)

  REAL*8,ALLOCATABLE:: pSLN(:,:)

  INTEGER,ALLOCATABLE:: pPLT(:,:)
  INTEGER,ALLOCATABLE:: pGRN(:,:)
  LOGICAL,ALLOCATABLE:: pEMR(:,:)

  REAL*8,ALLOCATABLE:: paGDH(:,:)
  REAL*8,ALLOCATABLE:: paVD(:,:)
  REAL*8,ALLOCATABLE:: pDVS(:,:)  !OZN
  REAL*8,ALLOCATABLE:: pDVSL(:,:)  !OZN
  REAL*8,ALLOCATABLE:: paO3FLX(:,:)  !OZN

  REAL*8,ALLOCATABLE:: pCDI(:,:)
  REAL*8,ALLOCATABLE:: pHDI(:,:)
  REAL*8,ALLOCATABLE:: pTAVE(:,:)
  REAL*8,ALLOCATABLE:: pNHED(:,:)


  REAL*8,ALLOCATABLE:: pLAI(:,:)
  REAL*8,ALLOCATABLE:: LAImx(:,:,:)
  REAL*8,ALLOCATABLE:: pDLFX(:,:)

  REAL*8,ALLOCATABLE:: pROT(:,:)
  REAL*8,ALLOCATABLE:: pHGT(:,:)

  REAL*8,ALLOCATABLE:: pWSH(:,:)
  REAL*8,ALLOCATABLE:: pWSO(:,:)
  REAL*8,ALLOCATABLE:: pWST(:,:)
  REAL*8,ALLOCATABLE:: pWLF(:,:)
  REAL*8,ALLOCATABLE:: pWRT(:,:)
  REAL*8,ALLOCATABLE:: pWAR(:,:)
  REAL*8,ALLOCATABLE:: pWIR(:,:)
  REAL*8,ALLOCATABLE:: pWDL(:,:)
  REAL*8,ALLOCATABLE:: pWGR(:,:)
  REAL*8,ALLOCATABLE:: pSWGR(:,:)
  REAL*8,ALLOCATABLE:: pNGR(:,:)
  REAL*8,ALLOCATABLE:: pNSP(:,:)

  REAL*8,ALLOCATABLE:: pGLF(:,:)
  REAL*8,ALLOCATABLE:: pGST(:,:)
  REAL*8,ALLOCATABLE:: pGRT(:,:)
  REAL*8,ALLOCATABLE:: pGSO(:,:)


  REAL*8,ALLOCATABLE:: TAVE5DAY(:) ! mean temperature during previous 5 days
  REAL*8,ALLOCATABLE:: pTAVE5DAY(:,:,:)
  INTEGER TAVE5CNT
  INTEGER,ALLOCATABLE:: pTAVE5CNT(:,:) 

!  [OUTPUT]
  REAL*8,ALLOCATABLE:: YLD(:,:,:)   ! Yield    [kg/ha]
  REAL*8,ALLOCATABLE:: PRM(:,:,:)   ! PRM - HI [-(%)]
  REAL*8,ALLOCATABLE:: OCDI(:,:,:)   ! CDI - HI [-(%)]

!  [FUNCTION]
  REAL*8 SINB
  REAL*8 CALHOUR


  
  !print *,"START"



  
  !! !!!!!!!!!!!!!!!!!!!! !!
  !!  <READ SETTING FILE> !!
  !! !!!!!!!!!!!!!!!!!!!! !!
  CALl RDSET(MODE,POINT,STYR,ENYR,STDOY,ENDOY,POINTLON,POINTLAT,TRES,clmTRES,WEST,EAST,NORTH,SOUTH,WERES,NSRES,IRR,&
             CRP_NAME,CRP_FILE,CNFILE,PLT_FILE,YLD_FILE,PRM_FILE,CDI_FILE,CLIM_HEAD,CLIM_FOOT,SMPL_FILE)

  NLON=1.D0
  NLAT=1.D0
  
!!!< Allocate >
  ALLOCATE(INTMP(NLON,NLAT,365*86400/clmTRES),INTMX(NLON,NLAT,365*86400/clmTRES),INTMN(NLON,NLAT,365*86400/clmTRES),INPRC(NLON,NLAT,365*86400/clmTRES),INRSD(NLON,NLAT,365*86400/clmTRES),INSHM(NLON,NLAT,365*86400/clmTRES),INWND(NLON,NLAT,365*86400/clmTRES),INPRS(NLON,NLAT,365*86400/clmTRES))
  ALLOCATE(INOZN(NLON,NLAT,365*86400/clmTRES))  !OZN

  ALLOCATE(PLTDOY(NLON,NLAT),GDHm(NLON,NLAT),GDHh(NLON,NLAT),NFERT(NLON,NLAT),Vsat(NLON,NLAT),LND(NLON,NLAT))
  ALLOCATE(SOILTXT(NLON,NLAT))
  ALLOCATE(BUF(NLON,NLAT))
  
  ALLOCATE(pPLT(NLON,NLAT),pGRN(NLON,NLAT),pEMR(NLON,NLAT))
  ALLOCATE(paGDH(NLON,NLAT))
  ALLOCATE(paVD(NLON,NLAT))
  ALLOCATE(pDVS(NLON,NLAT)) !OZN
  ALLOCATE(pDVSL(NLON,NLAT)) !OZN
  ALLOCATE(paO3FLX(NLON,NLAT))  !OZN
  ALLOCATE(pWSL(NLON,NLAT,NSL),pWSTRS(NLON,NLAT))
  ALLOCATE(pSLN(NLON,NLAT))

  ALLOCATE(pCDI(NLON,NLAT))
  ALLOCATE(pHDI(NLON,NLAT))
  ALLOCATE(pTAVE(NLON,NLAT))
  ALLOCATE(pNHED(NLON,NLAT))

  ALLOCATE(TAVE5DAY(86400/TRES*5))
  ALLOCATE(pTAVE5DAY(NLON,NLAT,86400/TRES*5)) !5days
  ALLOCATE(pTAVE5CNT(NLON,NLAT))

  
  ALLOCATE(pLAI(NLON,NLAT),pROT(NLON,NLAT),pHGT(NLON,NLAT))
  ALLOCATE(pDLFX(NLON,NLAT))
  ALLOCATE(pWSH(NLON,NLAT),pWSO(NLON,NLAT),pWST(NLON,NLAT),pWLF(NLON,NLAT),pWRT(NLON,NLAT),pWAR(NLON,NLAT),pWIR(NLON,NLAT),pWDL(NLON,NLAT),pGLF(NLON,NLAT),pGST(NLON,NLAT),pGRT(NLON,NLAT),pGSO(NLON,NLAT))
  ALLOCATE(pWGR(NLON,NLAT),pSWGR(NLON,NLAT),pNGR(NLON,NLAT),pNSP(NLON,NLAT))

  ALLOCATE(YLD(NLON,NLAT,ENYR-STYR+1),PRM(NLON,NLAT,ENYR-STYR+1),OCDI(NLON,NLAT,ENYR-STYR+1),LAImx(NLON,NLAT,ENYR-STYR+1))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! [CROP PARAMETER FILE] !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL RDPRM(CRP_FILE,RESPCP,EFFCON,ATHETA,BTHETA,MH2O,BH2O,KN,LTCH,ZKCA,ZKCB,ZKOA,ZKOB,GMMA,GMMB,RLFV,TLFV,RLFN,TLFN,hDVS,CFLF,CFST,CFRT,CFSO,VN,TB,TO,TH,LEFY0,LEFX1,LEFY1,LEFX2,LEFY2,LEFX3,LEFY3,PNCLX1,PNCLY1,PNCLX2,PNCLY2,PNCLX3,PNCLY3,DLFX1,DLFY1,DLFX2,DLFY2,DLFX3,DLFY3,LLFst,kLLF,RTX,RTY,RTX2,FSTR,SLWYA,SLWYB,SLWX,HGTAA,HGTAB,HGTBA,HGTBB,GZRT,MXRT,GMMSL,SLNX1,SLNX2,SLNX3,SLNYMX,SLNYMN,SLNK,TCmin,THcrit,HI,PLTDIF,HVT_TAVE)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! [SOIL TEXTURE, PLTDOY, and GDHm] !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(MODE=="ONE")THEN
      SOILTXT(NLON,NLAT)=7
      LND(:,:)=1.0
      WND_HGT=10.0
  END IF


  IF(CRP_NAME .eq. "Wheat")THEN

     Vsat=0.0

  END IF

  !print *,"SOIL, PLT, GDD OK"
  

  
  !! !!!!!!!!!!!!!!!!!!!!!!!!!! !!
  !! Initialize state variables !!
  !! !!!!!!!!!!!!!!!!!!!!!!!!!! !!
  pPLT(:,:)=0
  pGRN(:,:)=0
  pEMR(:,:)=.FALSE.

  paGDH(:,:)=0.D0
  paVD(:,:)=0.D0
  pDVS(:,:)=0.D0   !OZN
  pDVSL(:,:)=0.D0   !OZN
  paO3FLX(:,:)=0.D0 !OZN
  pLAI(:,:)=0.D0
  pROT(:,:)=0.D0
  pHGT(:,:)=0.D0
  pDLFX(:,:)=0.D0

  pCDI(:,:)=0.D0
  pHDI(:,:)=0.D0
  pTAVE(:,:)=0.D0
  pNHED(:,:)=0.D0

  pTAVE5DAY(:,:,:)=0.D0
  pTAVE5CNT(:,:)=0

  
  pWSL(:,:,:)=1.D0
  pWSTRS(:,:)=1.D0
  pSLN(:,:)=0.D0

  pWSH(:,:)=0.D0
  pWSO(:,:)=0.D0
  pWST(:,:)=0.D0
  pWLF(:,:)=0.D0
  pWRT(:,:)=0.D0
  pWAR(:,:)=0.D0
  pWIR(:,:)=0.D0
  pWDL(:,:)=0.D0
  pWGR(:,:)=0.D0
  pSWGR(:,:)=0.D0
  pNGR(:,:)=0.D0
  pNSP(:,:)=0.D0

  pGLF(:,:)=0.D0
  pGST(:,:)=0.D0
  pGRT(:,:)=0.D0
  pGSO(:,:)=0.D0
  
  !!!!! READING SMPL FILE !!!!!
  OPEN(20,file=SMPL_FILE,status='old',iostat=IER)
  !  print *,"koko",IER
  IF(IER .NE. 0)THEN
     print *, "SMPL.txt could not be opened"
     print *,SMPL_FILE
     STOP
  ELSE
   READ(20,*) SLNYMN
   READ(20,*) HI
   READ(20,*) TCmin
   READ(20,*) THcrit
   !READ(20,*) SLNY1
   !READ(20,*) GDHh(NLON,NLAT)
   READ(20,*) hDVS
   READ(20,*) GDHm(NLON,NLAT)
END IF
CLOSE(20)

SLNY1=1.D0

  !! !!!!!!!!!!!!!!!!!!! !!
  !! [Calculation START] !!
  !! !!!!!!!!!!!!!!!!!!! !!
  
  DO YEAR = STYR,ENYR
   !print *,YEAR
   
   !! !!!!!!!!!!!!!!!!!!!!! !!
   !! [READING YEARLY FILE] !!
   !! !!!!!!!!!!!!!!!!!!!!! !!
   
   IF(MODE=="ONE")THEN
      CALL RDCLIM(INTMN,CLIM_HEAD,POINT,'tmin ',YEAR,CLIM_FOOT,NLON,NLAT,365)
      !print *,"TMN OK"
      CALL RDCLIM(INTMX,CLIM_HEAD,POINT,'tmax ',YEAR,CLIM_FOOT,NLON,NLAT,365)
      !print *,"TMX OK"
      CALL RDCLIM(INPRC,CLIM_HEAD,POINT,'prc  ',YEAR,CLIM_FOOT,NLON,NLAT,365)
      !print *,"PRC OK"
      CALL RDCLIM(INRSD,CLIM_HEAD,POINT,'rsd  ',YEAR,CLIM_FOOT,NLON,NLAT,365)
      !print *,"RSD OK"
      CALL RDCLIM(INSHM,CLIM_HEAD,POINT,'shm  ',YEAR,CLIM_FOOT,NLON,NLAT,365)
      !print *,"SHM OK"
      CALL RDCLIM(INWND,CLIM_HEAD,POINT,'wnd  ',YEAR,CLIM_FOOT,NLON,NLAT,365)
      !print *,"WND OK"
      CALL RDCLIM(INPRS,CLIM_HEAD,POINT,'prs  ',YEAR,CLIM_FOOT,NLON,NLAT,365)
      !print *,"PRS OK"
      INOZN=0.0D0
      
      CALL RDCO2_NFERT(CO2PPM,NFERT(NLON,NLAT),CNFILE,YEAR)
      
   END IF
   
   !! <NITROGEN>
   
   DO I=1,NLON
      DO J=1,NLAT
         IF(NFERT(I,J) > 10000000.D0)THEN
            NFERT(I,J)=-1.D0
         END IF
         IF(LND(I,J) > 10000000.D0)THEN
            LND(I,J)=-1.D0
         END IF
      END DO
   END DO
   
   
   !! !!!!!!!!!!!!!! !!
   !! INITIALIZATION !!
   !! !!!!!!!!!!!!!! !!
   
   
   !! [Yield initialization] !!
     YLD(:,:,YEAR-STYR+1)=0.D0
     LAImx(:,:,YEAR-STYR+1)=0.D0
     
     CALL RDPLT(PLTDOY(NLON,NLAT),YEAR,POINT,PLT_FILE)
     
    PNCLY2 = 0.7D0 + 0.3D0 / (1960-1896) * (YEAR-1896)
    PNCLY2 = min(PNCLY2,1.D0)
    PNCLY3 = PNCLY2

     DO ILON=1,NLON


       !print *,ILON,"/",NLON

        DO ILAT=1,NLAT

          IF(MODE=="ONE")THEN
            LON = POINTLON
            LAT = POINTLAT
          ELSE
            LON = WEST + WERES * ILON - (WERES/2.D0)
            LAT = SOUTH + NSRES * ILAT - (NSRES/2.D0)
          END IF

!         print *,ILAT,LAT
           
         IF(LON > 180.D0)THEN
            LON = LON - 360.D0
         END IF
           
!                    print *,GDHm(ILON,ILAT,1,GRW),PLTDOY(ILON,ILAT,1,GRW),LON,LAT
                       
                    IF(LND(ILON,ILAT) > 0.D0 .AND. NFERT(ILON,ILAT) > -1.D0 .AND. &
                       PLTDOY(ILON,ILAT) < 1000 .AND. PLTDOY(ILON,ILAT) > 0 .AND. &
                       INPRS(ILON,ILAT,1) < 100000000.0 .AND. INPRS(ILON,ILAT,1) > 0.0 .AND. &
                       GDHm(ILON,ILAT) > 0.D0 .AND. SOILTXT(ILON,ILAT)<12 .AND. SOILTXT(ILON,ILAT) > 0)THEN !YM

                 

                       !! Update state variables !!!
                       PLT=pPLT(ILON,ILAT)
                       GRN=pGRN(ILON,ILAT)
                       EMR=pEMR(ILON,ILAT)
                       aGDH=paGDH(ILON,ILAT)
                       aVD=paVD(ILON,ILAT)
                       DVS=pDVS(ILON,ILAT)   !OZN
                       DVSL=pDVSL(ILON,ILAT)   !OZN
                       aO3FLX=paO3FLX(ILON,ILAT)  !OZN
                       LAI=pLAI(ILON,ILAT)
                       ROT=pROT(ILON,ILAT)
                       HGT=pHGT(ILON,ILAT)
                       DLFX=pDLFX(ILON,ILAT)
                       

                       CDI=pCDI(ILON,ILAT)
                       HDI=pHDI(ILON,ILAT)
                       TAVE=pTAVE(ILON,ILAT)
                       
                       NHED=pNHED(ILON,ILAT)
                       
                       TAVE5DAY(:)=pTAVE5DAY(ILON,ILAT,:)
                       TAVE5CNT=pTAVE5CNT(ILON,ILAT)
                       
                       
                       WSH=pWSH(ILON,ILAT)
                       WSO=pWSO(ILON,ILAT)
                       WST=pWST(ILON,ILAT)
                       WLF=pWLF(ILON,ILAT)
                       WRT=pWRT(ILON,ILAT)
                       WAR=pWAR(ILON,ILAT)
                       WIR=pWIR(ILON,ILAT)
                       WDL=pWDL(ILON,ILAT)
                       WGR=pWGR(ILON,ILAT)
                       SWGR=pSWGR(ILON,ILAT)
                       NGR=pNGR(ILON,ILAT)
                       NSP=pNSP(ILON,ILAT)
                       
                       GLF=pGLF(ILON,ILAT)
                       GST=pGST(ILON,ILAT)
                       GRT=pGRT(ILON,ILAT)
                       GSO=pGSO(ILON,ILAT)
                       
                       
                       SLN=pSLN(ILON,ILAT)
                       WSL(:)=pWSL(ILON,ILAT,:)
                       WSTRS=pWSTRS(ILON,ILAT)
          
                       IF(YEAR==STYR .AND. YEAR==ENYR)THEN
                          STDD=STDOY
                          ENDD=ENDOY
                       ELSE if(YEAR==STYR)THEN
                          STDD=STDOY
                          ENDD=365
                       ELSE if(YEAR==ENYR)THEN
                          STDD=1
                          ENDD=ENDOY
                       ELSE
                          STDD=1
                          ENDD=365
                       END IF
                       
                          

                       DO DOY  = STDD,ENDD 
                        !print *,DOY
                        !print *,HGT 


                       TMX=INTMX(ILON,ILAT,DOY)
                       
                       INTSINB = 0.D0
                       DO IHOUR=1,86400/TRES
                          HOUR = CALHOUR(IHOUR,TRES)
                          INTSINB = INTSINB + SINB(DOY,HOUR,LAT)*DBLE(TRES)
                       END DO
                       
                       PDOY = DOY - 1 
                       NDOY = DOY + 1
                       IF(PDOY==0)THEN
                          PDOY=365
                       END IF
                       IF(NDOY==366)THEN
                          NDOY=1
                       END IF
                       
                       DO IHOUR=1,86400/TRES
                          
                          HOUR =  CALHOUR(IHOUR,TRES)    
                          !print *,HOUR

                          IF(clmTRES<86400)THEN
                             TMP=INTMP(ILON,ILAT,(DOY-STDD)*(86400/TRES)+IHOUR)
                             PRC=INPRC(ILON,ILAT,(DOY-STDD)*(86400/TRES)+IHOUR)
                             RSD=INRSD(ILON,ILAT,(DOY-STDD)*(86400/TRES)+IHOUR)
                             SHM=INSHM(ILON,ILAT,(DOY-STDD)*(86400/TRES)+IHOUR)
                             WND=INWND(ILON,ILAT,(DOY-STDD)*(86400/TRES)+IHOUR)
                             PRS=INPRS(ILON,ILAT,(DOY-STDD)*(86400/TRES)+IHOUR)
                             OZN=INOZN(ILON,ILAT,(DOY-STDD)*(86400/TRES)+IHOUR)
                          ELSE
                             
                             CALL TINTERP(TMP,PRC,RSD,SHM,WND,PRS,OZN,INTMX(ILON,ILAT,PDOY)-273.15D0,INTMX(ILON,ILAT,DOY)-273.15D0,INTMX(ILON,ILAT,NDOY)-273.15D0,INTMN(ILON,ILAT,PDOY)-273.15D0,INTMN(ILON,ILAT,DOY)-273.15D0,INTMN(ILON,ILAT,NDOY)-273.15D0,INPRC(ILON,ILAT,DOY),INRSD(ILON,ILAT,DOY),INSHM(ILON,ILAT,DOY),INWND(ILON,ILAT,DOY),INPRS(ILON,ILAT,DOY),INOZN(ILON,ILAT,DOY),PDOY,DOY,NDOY,HOUR,LAT,TRES,INTSINB,WND_HGT)
                             
                          END IF
                          
                          WND=MAX(WND,0.001D0)

                          IF(TAVE5CNT < (86400 / TRES * 5))THEN
                             TAVE5CNT = TAVE5CNT + 1
                             TAVE5DAY(TAVE5CNT) = TMP - 273.15D0
                          ELSE
                             TAVE5DAY(1:(TAVE5CNT-1))=TAVE5DAY(2:TAVE5CNT)
                             TAVE5DAY(TAVE5CNT)=TMP - 273.15d0
                          END IF
                          
                          
                          !! Absorbed PAR, Vmax, in Sunlit and Shade leaves per LAI
                          CALL RAD(QPARSNLF,QPARSHLF,VMXSNLF,VMXSHLF,LAISN,LAISH,SLN,KN,RSD,LAI,RLFV,TLFV,RLFN,TLFN,LAT,DOY,HOUR,ILON,ILAT,CRP_NAME)
                          !print *,"LAI",LAI
                          !!print *,"LATSN",LAISN
                          !print *,"LAISH",LAISH

                          !! Canopy Gross Assimilation and Respiration, Bulk coefficient of conductance for vapor
                          CALL PHSYN(GPP,RSP,TSP,QPARSNLF,QPARSHLF,VMXSNLF,VMXSHLF,LAISN,LAISH,TMP,SHM,PRS,WND,CO2PPM,OZN,WSTRS,aO3FLX,LAI,HGT,DVS,hDVS,RESPCP,EFFCON,ATHETA,BTHETA,MH2O,BH2O,LTCH,ZKCA,ZKCB,ZKOA,ZKOB,GMMA,GMMB,TRES) 
                          
                          !! Crop simulation
                          CALL CROP(YLD(ILON,ILAT,YEAR-STYR+1),PRM(ILON,ILAT,YEAR-STYR+1),OCDI(NLON,NLAT,YEAR-STYR+1),PLT,EMR,GRN,DVS,DVSL,aGDH,aVD,CDI,HDI,TAVE,TMX,NHED,aO3FLX,LAI,LAImx(ILON,ILAT,YEAR-STYR+1),HGT,ROT,WSH,WSO,WST,WLF,WRT,WAR,WIR,WDL,WGR,SWGR,NGR,NSP,GLF,GST,GRT,GSO,SLN,NFERT(ILON,ILAT),Vsat(ILON,ILAT),GPP,RSP,TMP,PLTDOY(ILON,ILAT)+INT(PLTDIF),DOY,HOUR,TRES,GDHm(ILON,ILAT),hDVS,CFLF,CFST,CFRT,CFSO,VN,TB,TO,TH,LEFY0,LEFX1,LEFY1,LEFX2,LEFY2,LEFX3,LEFY3,PNCLX1,PNCLY1,PNCLX2,PNCLY2,PNCLX3,PNCLY3,DLFX1,DLFY1,DLFX2,DLFY2,DLFX3,DLFY3,RTX,RTY,RTX2,FSTR,SLWYA,SLWYB,SLWX,HGTAA,HGTAB,HGTBA,HGTBB,GZRT,MXRT,SLNX1,SLNY1,SLNX2,SLNX3,SLNYMX,SLNYMN,SLNK,WSTRS,ILON,ILAT,IRR,TCmin,THcrit,HI,TAVE5DAY,TAVE5CNT,HVT_TAVE,LLFst,kLLF,DLFX,CO2PPM,YEAR) 

                          !! Soil water balance
                          CALL SOIL(WSTRS,WSL,TSP*EL,PRC,ROT,IRR,GMMSL,TRES,SOILTXT(ILON,ILAT),TMP,PRS,WND,SHM,HGT,PLT,LON,LAT,CRP_NAME) 

                        END DO

                        !!print *,"DAYEND",HGT
                        !print *,"DAYEND",LAI
                       
                    END DO
                       
                       !! Store state variables !!!!!!!!
                       pPLT(ILON,ILAT)=PLT
                       pGRN(ILON,ILAT)=GRN
                       pEMR(ILON,ILAT)=EMR
                       
                       paGDH(ILON,ILAT)=aGDH
                       paVD(ILON,ILAT)=aVD

                       pDVS(ILON,ILAT)=DVS   !OZN
                       pDVSL(ILON,ILAT)=DVSL   !OZN
                       paO3FLX(ILON,ILAT)=aO3FLX !OZN
                       pLAI(ILON,ILAT)=LAI
                       pROT(ILON,ILAT)=ROT
                       pHGT(ILON,ILAT)=HGT
                       pDLFX(ILON,ILAT)=DLFX
                       
                       pCDI(ILON,ILAT)=CDI
                       pHDI(ILON,ILAT)=HDI
                       pTAVE(ILON,ILAT)=TAVE
                       pNHED(ILON,ILAT)=NHED
                       pTAVE5DAY(ILON,ILAT,:)=TAVE5DAY(:)
                       pTAVE5CNT(ILON,ILAT)=TAVE5CNT
                       
                       
                       pWSH(ILON,ILAT)=WSH
                       pWSO(ILON,ILAT)=WSO
                       pWST(ILON,ILAT)=WST
                       pWLF(ILON,ILAT)=WLF
                       pWRT(ILON,ILAT)=WRT
                       pWAR(ILON,ILAT)=WAR
                       pWIR(ILON,ILAT)=WIR
                       pWDL(ILON,ILAT)=WDL
                       pWGR(ILON,ILAT)=WGR
                       pSWGR(ILON,ILAT)=SWGR
                       pNGR(ILON,ILAT)=NGR
                       pNSP(ILON,ILAT)=NSP
                       
                       pGLF(ILON,ILAT)=GLF
                       pGST(ILON,ILAT)=GST
                       pGRT(ILON,ILAT)=GRT
                       pGSO(ILON,ILAT)=GSO
                       
                       
                       pSLN(ILON,ILAT)=SLN
                       pWSL(ILON,ILAT,:)=WSL(:)
                       pWSTRS(ILON,ILAT)=WSTRS
                       !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    END IF
              
                 END DO !ILAT
              END DO !ILON

              ! tentative yield
              CALL WRTXT_YLD(STYR,ENYR,YLD_FILE,YLD(NLON,NLAT,:),NLON,NLAT) !YLD
              CALL WRTXT_YLD(STYR,ENYR,PRM_FILE,PRM(NLON,NLAT,:),NLON,NLAT) !HI
              CALL WRTXT_YLD(STYR,ENYR,CDI_FILE,OCDI(NLON,NLAT,:),NLON,NLAT)
              
            END DO            ! YEAR
         !print *,YLD(NLON,NLAT,:)
         !print *,PRM(NLON,NLAT,:)

END PROGRAM MATCRO
