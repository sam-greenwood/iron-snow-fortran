  module parameters_core

    implicit none

    integer         , parameter :: Nc=1000, Ni=1000, Ns=1000, N=Nc+Ni+Ns 
    double precision, parameter :: hl=1.248e3          !Half life of K40
    double precision, parameter :: dmudT0=-0.00023, dmudT1=-6.32184e-08 

    double precision, parameter :: rho_zeroP = 6.0e3   !Density, zero pressure
    double precision, parameter :: K_zeroP   = 500e9   !Bulk modulus, zero pressure

    integer, parameter          :: AO         = 16    , AS         = 32    , ASi         = 28,   AFe=56
    double precision, parameter :: lambdaO_oc = 3.25d0, lambdaS_oc = 6.15d0, lambdaSi_oc = 3.6d0
    double precision, parameter :: lambdaO_ic = 0.00d0, lambdaS_ic = 5.9d0 , lambdaSi_ic = 2.7d0
    double precision, parameter :: dmuO       =-2.6d0 , dmuS       =-0.25d0, dmuSi       =-0.05d0
    double precision, parameter :: DO         = 1e-8  , DS         = 5e-9  , DSi         = 5e-9
    double precision, parameter :: alphac_cO  = 1.1d0 , alphac_cS  = 0.64  , alphac_cSi  = 0.87, alphaT_c=585d-7
    double precision, parameter :: Rh=-27.7d6

  contains

  end module parameters_core
