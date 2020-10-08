  module parameters_core

    implicit none

    integer         , parameter :: Nc=400, Ni=200, N=Nc+Ni  !Number of pts in the outer, inner core & total
    double precision, parameter :: rc=3480e3     !CMB radius
    double precision, parameter :: hl=1.248e3    !Half life of K40
    double precision, parameter :: dmudT0=-0.00023, dmudT1=-6.32184e-08 

    double precision, parameter :: rho_zeroP = 7.9e3 !Density, zero pressure
    double precision, parameter :: K_zeroP   = 500e9   !Bulk modulus, zero pressure

    integer, parameter          :: AFe=56, AO=16, AS=30, ASi=28
    double precision, parameter :: lambdaO_oc = 3.25d0, lambdaS_oc = 6.2d0, lambdaSi_oc = 3.5d0
    double precision, parameter :: lambdaO_ic = 0.00d0, lambdaS_ic = 5.9d0, lambdaSi_ic = 0.0d0
    double precision, parameter :: dmuO=-2.6d0, dmuS=-0.25d0
    double precision, parameter :: DO         = 1e-8  , DS         = 5e-9 , DSi         = 5e-9
    double precision, parameter :: alphac_cS=0.64, alphac_cSi=0.87
    double precision, parameter :: Rh=-27.7d6

  contains

  end module parameters_core

