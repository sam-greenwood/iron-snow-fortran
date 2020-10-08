  module mantle_evolution

    use parameters_mantle

    implicit none

    double precision :: Tm
    double precision :: rm
    double precision :: cp_m 
    double precision :: rho_m
    double precision :: g_m 
    double precision :: k_m  
    double precision :: kappa_m
    double precision :: alpha_m
    double precision :: Ts   
    double precision :: gamma
    double precision :: eta0 
    double precision :: A_m
    double precision :: afac, bfac, Rac
    character (len=100) :: out_file
    integer             :: out_file_len

  contains

  subroutine read_input_mantle(fname)

    character (len=200)   :: line, fname

    open(5,file=fname)

101   continue
      read(5, 80) line
      if(line(1:1) .eq. '*') goto 101
      out_file     = adjustl(line)
      out_file_len = index(out_file, ' ')-1
      out_file    = out_file(1:out_file_len)
 
102   continue
      read(5, 80) line
      if(line(1:1) .eq. '*') goto 102
      read(line, *) rho_m, g_m, rm
      A_m     = 4d0*pi*rm**2

103   read(5, 80) line
      if(line(1:1) .eq. '*') goto 103
      read(line, *) cp_m, k_m, alpha_m, eta0
      kappa_m = k_m/(rho_m*cp_m)
 
104   read(5, 80) line
      if(line(1:1) .eq. '*') goto 104
      read(line, *) Ts, Tm, gamma

105   read(5, 80) line
      if(line(1:1) .eq. '*') goto 105
      read(line, *) afac, bfac, Rac

 80   FORMAT(A)

  end subroutine read_input_mantle

  subroutine open_files_mantle(char)

    character (len=4)   :: char
    character (len=100) :: energy_file_tot, energy_file_pua 

    energy_file_tot  = out_file(1:out_file_len)//"_"//char//'_energy_tot'
    energy_file_pua  = out_file(1:out_file_len)//"_"//char//'_energy_pua'

    open(166,file=energy_file_tot)
    open(167,file=energy_file_pua)
    write(166,'(10A16)') 'Time (Myr)','Qm','Qc','Qr','dTmdt','Tm','deltau','deltab','etai'
    write(167,'(10A16)') 'Time (Myr)','Qm','Qc','Qr','dTmdt','Tm','deltau','deltab','etai'

  end subroutine open_files_mantle

  subroutine mantle_evol_calc(Tc,rc,time,forback,Tm,Qc)

    double precision, intent(out):: Tm, Qc
    double precision, intent(in) :: Tc, rc, time
    integer         , intent(in) :: forback

    double precision :: V_m, M_m, A_c
    double precision :: Ti, etai, dTl
    double precision :: Fsl, Fc, fac, deltau, deltab, Hint
    double precision :: t, dTmdt

    A_c  = 4d0*pi*rc**2                    !Core area
    V_m  = 4d0*pi*(rm**3-rc**3)/3d0        !Mantle volume
    M_m  = V_m*rho_m                       !Mantle mass (const rho)

    !CD - Not clear what is the "interior" temperature. 
    !     Mantle has 1 temperature Tm so assume Ti=Tm
    !     Also note that eta is calibrated to viscosity at 1300 C
    Ti   = Tm!(Tc - Ts)/2d0                !Interior temperature
    etai = eta0 * dexp(-gamma*(Ti-1573.0)) !Interior viscosity

    fac    = rho_m * g_m *alpha_m / (kappa_m * etai)
    Fsl    = k_m / afac * fac**(1d0/3d0) * gamma**(-4d0/3d0)  !NS eq 3
    deltau = (Rac * afac**4d0 / fac / (8d0/gamma))**(1d0/3d0) !NS eq 6
   
    dTl    = Tc-Tm
    if(dTl < 0d0) dTl = -1d0*dTl
    deltab = 0.5*deltau*(gamma*dTl)**(-1d0/3d0) * dexp(-gamma*(Tc-Tm)/6d0)
                                                  !NS eq 7
    Fc     = k_m * (Tc-Tm) / deltab                        !NS eq 5

    ! Radiogenic heating
    t    = (ttot/1e6) - time
    Hint = (ppmK40*2.08e-4*hK40*dexp(log(2d0)*t/hlK40) + & 
            ppmU  *0.9927 *hU  *dexp(log(2d0)*t/hlU)   + & 
            ppmT  *        hT  *dexp(log(2d0)*t/hlT)) * M_m

    dTmdt= (-A_m*Fsl + A_c*Fc + Hint)/(cp_m*M_m)

    Qc   = A_c*Fc

    write(166,'(9E16.8)') time, A_m*Fsl, A_c*Fc, Hint, dTmdt*secingyr, Tm, deltau, deltab, etai
    write(167,'(9E16.8)') time, Fsl*1e3, Fc*1e3, Hint/A_m*1e3, dTmdt*secingyr, Tm, deltau, deltab, etai

    if(forback==1) Tm  = Tm + dTmdt*secinyr*dt        !New mantle temperature
    if(forback==2) Tm  = Tm - dTmdt*secinyr*dt        !New mantle temperature

  end subroutine mantle_evol_calc

! Tm = mantle T at base of upper BL
! Tc = core temperature
! Tb = T at top of lower BL
! Tl = T at top of upper BL
! dTl = Tc-Tb; lower BL T drop
! dTu = Tl-Tm; upper BL T drop
  subroutine mantle_evol_calc_breuer(Tc,rc,time,forback,Tm,Qc)

    double precision, intent(out):: Tm, Qc
    double precision, intent(in) :: Tc, rc, time
    integer         , intent(in) :: forback

    double precision :: V_m, M_m, A_c, rl, A_sl
    double precision :: Tb, Tl, Tref, etau, etal, Tli
    double precision :: dTu, dTl, deltau, deltal
    double precision :: Fsl, Fc, Hint
    double precision :: Racl, Racu, Rau, Ral
    double precision :: t, dTmdt
    double precision :: Ag

    Ag   = gamma                           !Activation energy
    Tref = 1600.0                          !Reference visc
    rl   = rm-350d3                        !ASSUME r_lid doesn't vary in time
    Racu = Rac                             !Crit Ra for upper TBL

    A_sl = 4d0*pi*rl**2                    !Area of base of lid
    A_c  = 4d0*pi*rc**2                    !Area of CMB
    V_m  = 4d0*pi*(rl**3-rc**3)/3d0        !Mantle volume
    M_m  = V_m*rho_m                       !Mantle mass (const rho)

    Tl   = Tm - afac * Rg*Tm**2/Ag
    Tb   = Tm !+ alpha_m*Tm*g*dR/cp_m      !T drop adiabatic increase, see below14
    Tli  = (Tc+Tm)/2d0

    dTu  = Tm-Tl !Chk - their stuff below 14 cant be right!
    dTl  = Tc-Tb

    etau  = eta0 * dexp((Ag/Rg)*(1d0/Tm -1d0/Tref))
    etal  = eta0 * dexp((Ag/Rg)*(1d0/Tli-1d0/Tref))   

    ! Assumes g and alpha and kappa are constant
    Rau   = alpha_m*rho_m*g_m*dTu*(rl-rc)**3 / (kappa_m * etau)
    Ral   = alpha_m*rho_m*g_m*dTl*(rl-rc)**3 / (kappa_m * etal)
    if(Ral<0.0) Ral = -1d0 * Ral
    Racl = 0.28*(alpha_m*rho_m*g_m*(Tc-Ts)*(rm-rc)**3 / (kappa_m * etal))**0.21  
    ! Thiriet never say where eta is evaluated in Ra_int
    ! I also assume dT and D here refer to whole shell. 

    deltau = (rl-rc)*(Racu/Rau)**bfac
    deltal = (rl-rc)*(Racl/Ral)**bfac

print *, Racl, Ral, rl-rc

    Fsl    = k_m * dTu / deltau
    Fc     = k_m * dTl / deltal

    ! Radiogenic heating
    t    = (ttot/1e6) - time
    Hint = (1.9d0*5.56e-13*dexp(log(2d0)*t/hlK40) + &
            1.50e-12*dexp(log(2d0)*t/hlU)   + &
            6.46e-14*dexp(log(2d0)*t/hlU235)+ &
            0.875*1.69e-12*dexp(log(2d0)*t/hlT)) * M_m

    dTmdt= (-A_sl*Fsl + A_c*Fc + Hint)/(cp_m*M_m) !GB08 eq 1

    Qc   = A_c*Fc

    write(166,'(9E16.8)') time, A_sl*Fsl, A_c*Fc, Hint, dTmdt*secingyr, Tm,deltau, deltal, etal
    write(167,'(10E16.8)') time, Fsl*1e3, Fc*1e3, Hint/A_sl*1e3, dTmdt*secingyr, Tm, deltau, deltal, etau, Racl

    if(forback==1) Tm  = Tm + dTmdt*secinyr*dt        !New mantle temperature
    if(forback==2) Tm  = Tm - dTmdt*secinyr*dt        !New mantle temperature

  end subroutine mantle_evol_calc_breuer

  end module mantle_evolution
