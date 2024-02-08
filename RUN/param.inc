
REAL*8,PARAMETER::WMH2O= 0.018015D0
REAL*8,PARAMETER::WMCO2 = 44.D-3
REAL*8,PARAMETER::TQICE = 273.15D0
INTEGER,PARAMETER::TNMAX = 2*24*365*10 
INTEGER,PARAMETER::YRMAX = 3000
REAL*8,PARAMETER::PO2A=20900.D0
REAL*8,PARAMETER::SECPD=86400.D0
REAL*8,PARAMETER::CO2CHO=30.D0/44.D0
REAL*8,PARAMETER::CO2C=12.D0/44.D0
REAL*8,PARAMETER::ST2GLU=1.111D0
REAL*8,PARAMETER::LSSGLU=0.947D0
REAL*8,PARAMETER::C2CHO=2.5D0
REAL*8,PARAMETER::CHO2C=0.4D0

REAL*8,PARAMETER::CFSR2GL=1.11D0


REAL*8,PARAMETER::CPGLF=0.408D0
REAL*8,PARAMETER::CPGST=0.365D0
REAL*8,PARAMETER::CPGRT=0.365D0
REAL*8,PARAMETER::CPGSO=0.357D0

INTEGER,PARAMETER::NSL=5        ! Number of soil layer
REAL*8,PARAMETER:: ZSL(NSL) = (/ 0.05D0, 0.2D0, 0.75D0 , 1.D0, 2.D0    /) !Thickness of each soil layer

INTEGER,PARAMETER::NCMAX=200
REAL*8,PARAMETER::DNCNPY=0.1D0

REAL*8,PARAMETER::ZA=2.D0
REAL*8,PARAMETER::ALBS=0.1D0
REAL*8,PARAMETER::EMIS=0.96D0



!! Physical parameter
      REAL*8     ER                          !! radius of the earth
      REAL*8     GRAV                        !! acceleration of gravity
      REAL*8     CP                          !! specific heat of air
      REAL*8     RAIR                        !! gas constant of air
      REAL*8     EL                          !! latent heat of cond.
      REAL*8     CPVAP                       !! specific heat of vapor
      REAL*8     RVAP                        !! gas constant of q
      REAL*8     DWATR                       !! density of water
      REAL*8     EPSV                        !! molec. weight water/air
      REAL*8     ES0                         !! saturate e 0deg C(Pa)
      REAL*8     STB                         !! Stefan-Boltzman
      REAL*8     FKARM                       !! Karman constant
      REAL*8     EPSVT                       !! 1/epsilon-1
      REAL*8     EMELT                       !! latent heat of melt.
      REAL*8     TMELT                       !! melting point of water
      REAL*8     CPWATR                      !! specific heat of water
      REAL*8     CPICE                       !! specific heat of ice
      REAL*8     TFRZS                       !! freeze point of sea

      REAL*8     KTW     !! thermal conductivity of water
      REAL*8     CPMN      !! specific heat of soil minerals

!
      PARAMETER ( ER     = 6370.D+3 )
      PARAMETER ( GRAV   = 9.8D0  )
      PARAMETER ( CP     = 1004.6D0  )
      PARAMETER ( RAIR   = 287.04D0  )
      PARAMETER ( EL     = 2.5D+6 )
      PARAMETER ( CPVAP  = 1810.D0  )
      PARAMETER ( RVAP   = 461.5D0  )
      PARAMETER ( DWATR  = 1000.D0  )
      PARAMETER ( ES0    = 611.D0  )
      PARAMETER ( STB    = 5.67D-8 )
      PARAMETER ( FKARM  = 0.4D0  )
      PARAMETER ( EMELT  = 3.4D+5 )
      PARAMETER ( TMELT  = 273.15D0 )
      PARAMETER ( CPWATR = 4200.D0 )
      PARAMETER ( CPICE  = 2000.D0 )
      PARAMETER ( TFRZS  = 271.35D0 )
      PARAMETER ( EPSV   = RAIR / RVAP )
      PARAMETER ( EPSVT  = 1.0D0/EPSV - 1.0D0 )
      PARAMETER(KTW = 0.6D0)
      PARAMETER(CPMN=870.D0)   !J/Kg/K


      REAL*8 HOURTMX   !! Hour at maximum temp
      REAL*8 NocT      !! Nocturnal time coefficient
      PARAMETER(HOURTMX=1.5D0)
      PARAMETER(NocT=4.D0)




!! mathematical 
   REAL*8 PI
   PARAMETER(PI=ACOS(-1.D0))


!! vege param
   REAL*8,PARAMETER::LFOR=0.5D0
   REAL*8,PARAMETER::STMIN=1.D-3


  REAL*8,PARAMETER::WCMAXR=0.1D-3