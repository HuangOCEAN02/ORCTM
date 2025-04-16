module chpar
implicit none
! type setting
INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6,37)
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)
INTEGER, PARAMETER :: wp = dp   ! working precision
INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)
INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)
! dimension of regional grid
integer,parameter::IE=817 
integer,parameter::JE=619
integer,parameter::KE=40
! dimension of total Grid numbers
integer,parameter::ITO=IE*2+1
integer,parameter::JTO=JE*2+1
! first level
real,parameter::lev1st=3.5702  ! 3.5702 from BRAVE60
! max depth (fill the trenches)
real,parameter::max_dep=6100
! dimension of climate TS data WOA18
integer,parameter::IDE=1440,JDE=720,KDE=57
! dimension of FVCOM INITIAL FIELDS 
integer,parameter::IDE_G=561,JDE_G=521,KIE=16
! resolution of climate TS data WOA18
real,parameter::RESOL=0.25
! dimension of JRA-55-TL319 forcing data
integer,parameter::IFE=640,JFE=320
! resolution of JRA-55-TL319 forcing data
real,parameter::resol_fe=0.5625
! input days setting
integer,parameter::LDAY=365,LMON=12,LDAY_NO_CLI=1
! restart setting
integer,parameter::FRONT_END=0
! FORCE_DWLW
integer,parameter::FORCE_DWLW=1
! StateEquation 
integer,parameter :: whichEOS=1
! others

! earth parameters
REAL(dp), PARAMETER :: api   = 3.1415926
REAL(dp), PARAMETER :: aradtogra = 180./api

!>  RADIUS OF THE ABOVE MENTIONED PLANET
REAL(wp), PARAMETER :: radius = 6371.E+3_wp
!> time per rotation in seconds, sidereal day
REAL(wp), PARAMETER :: trot = 86164._wp
!>  ANGULAR VELOCITY OF OUR NICE BLUE PLANET
REAL(wp), PARAMETER :: OMEGA = 2._wp * api/trot

end module chpar
