  module core_model

    use parameters
    use parameters_core

    implicit none

    character (len=100) :: out_file
    double precision :: rho_cen, Pc, cp_c
    double precision :: k0_oc, k1_oc, k2_oc
    double precision :: Tm0_c, Tm1_c, Tm2_c, Tm3_c, Tm4_c, meltdep
    double precision :: M_c, M_oc, M_ic, M_s
    double precision :: rho0_oc, rho1_oc, rho2_oc, rho3_oc, rho0_ic, rho2_ic
    double precision :: T_cen,  t1_c,  t2_c,  t3_c
    double precision :: ds0_c, ds1_c, ds2_c, ds3_c, ds4_c
    double precision :: pr_c(0:N)   , rhor_c(0:N) , gr_c(0:N) , psir_c(0:N)
    double precision :: Tar_c(0:N)  , Tm_c(0:N)   , ds_c(0:N)
    double precision :: dTadr_c(0:N), dTmdP_c(0:N), kr_c(0:N)
    double precision ::  r_c(0:N),r2_c(0:N),r3_c(0:N),r4_c(0:N),r5_c(0:N)
    double precision :: r6_c(0:N),r7_c(0:N),r8_c(0:N),r9_c(0:N),r10_c(0:N),r11_c(0:N)
    double precision :: ri0, rs, h0_c, h, rc
    double precision, allocatable :: cmbflux(:), time_file(:), dt_file(:)
    integer          :: nt, out_file_len, sol, ah, hor
    integer          :: iteration

  contains

!******************************************************************************
! read_input: reads the input file davies_TH_params
!
! Inputs: Q  - CMB heat flux
!         Tc - CMb temperature
!         ri - Inner core radius
!         EJ - Ohmic heating
!
! Comments in the input file are lines starting with an *
! If the solution method sol=0 then Q is read in from a file
!******************************************************************************
  subroutine read_input_core(fname,Q,Tc,ri,EJ,cbarO_oc,cbarS_oc,cbarSi_oc)

    double precision      :: cbarO_oc, cbarS_oc, cbarSi_oc
    double precision      :: tmp, Tc, ri, Q, EJ
    double precision      :: Qpres, Qfac, Qtau
    character (len=200)   :: line, fname

!   This was used to read in the filename containing input variables.
!   Now this name is hard-coded

    open(5,file=fname)

21    continue
      read(5, 80) line
      if(line(1:1) .eq. '*') goto 21
      out_file     = adjustl(line)
      out_file_len = index(out_file, ' ')-1
      out_file    = out_file(1:out_file_len)

22    continue
      read(5, 80) line
      if(line(1:1) .eq. '*') goto 22
      read(line, *) ah, hor

23    continue
      read(5, 80) line
      if(line(1:1) .eq. '*') goto 23
      read(line, *) rho_cen, Pc, cp_c

24    continue
      read(5, 80) line
      if(line(1:1) .eq. '*') goto 24
      read(line, *) rho0_oc, rho1_oc, rho2_oc, rho3_oc, rho0_ic, rho2_ic

25    continue
      read(5, 80) line
      if(line(1:1) .eq. '*') goto 25
      read(line, *) ds0_c, ds1_c, ds2_c, ds3_c, ds4_c

26    continue
      read(5, 80) line
      if(line(1:1) .eq. '*') goto 26
      read(line, *) Tm0_c, Tm1_c, Tm2_c, Tm3_c, Tm4_c, meltdep

27    continue
      read(5, 80) line
      if(line(1:1) .eq. '*') goto 27
      read(line, *) T_cen, t1_c, t2_c, t3_c

28    continue
      read(5, 80) line
      if(line(1:1) .eq. '*') goto 28
      read(line, *) k0_oc, k1_oc, k2_oc

29    continue
      read(5, 80) line
      if(line(1:1) .eq. '*') goto 29
      read(line, *) cbarO_oc, cbarS_oc, cbarSi_oc

30    continue
      read(5, 80) line
      if(line(1:1) .eq. '*') goto 30
      read(line, *) Tc, ri, rc, h0_c
      h = h0_c
      rs= rc

31    continue
      read(5, 80) line
      if(line(1:1) .eq. '*') goto 31
      read(line, *) dt,ttot

      Q  = 0.d0
      Ej = 0.d0
32    continue
      read(5, 80) line
      if(line(1:1) .eq. '*') goto 32
      read(line, *) sol, tmp

      if(sol==1 .or. sol==3) then
        Q  = tmp
      elseif(sol==2) then
        EJ = tmp
      endif

!     If read in Q, nt to be determined
      nt = 0
      if(sol==0) then

33      continue
        read(5, 80) line
        if(line(1:1) .eq. '*') goto 33

        call read_cmbhf_core(line)
      else
        nt = int(ttot/dt)

        allocate(cmbflux(0:nt))

        cmbflux = 0.d0
      endif

      ! Read in variables for functional form of Q
      if(sol==4) then 
34      continue
        read(5, 80) line
        if(line(1:1) .eq. '*') goto 34
        read(line, *) Qpres, Qfac, Qtau
        cmbflux(0) = Qpres; cmbflux(1) = Qfac; cmbflux(2) = Qtau
      endif

      write(*,*) '# timepts = ', nt

      close(5)

 80   FORMAT(A)

  end subroutine read_input_core

!******************************************************************************
! read_cmbhf: reads the CMB heat flux Q from a file
!
! Inputs: Qcmb_file - name of file containing Q
!                     format is 2 columns: t, Q
!
! n is a parameter set arbitrarily large (hopefully larger than the # rows!)
!******************************************************************************
  subroutine read_cmbhf_core(Qcmb_file)

    integer         , parameter   :: n=100000000
    double precision              :: tmp(n), t(n)
    character (len=200) :: Qcmb_file
    integer :: i

    open(99,file=Qcmb_file)

    t  = 0.d0; tmp = 0.d0
    nt = 0

    do i = 0, n
      read(99,*,end=100) t(i), tmp(i)
      nt = nt + 1
    enddo

100 continue

    nt = nt - 1

    allocate(cmbflux(0:nt))
    allocate(time_file(0:nt))
    allocate(dt_file(0:nt))

!   Set up time and Q so they run backwards.
    dt_file = 0.d0
    time_file = 0.d0
    cmbflux = 0.d0
    do i = 0, nt
      time_file(i) = t(nt-i)*1e3
      cmbflux(i)   = tmp(nt-i)*1e12                                  !Q in TW
      if(i > 0) dt_file(i) = dabs(time_file(i) - time_file(i-1))*1e6 !dt in Myr
!      write(546,*) time_file(i), dt_file(i-1), cmbflux(i)
    enddo

  end subroutine read_cmbhf_core

!******************************************************************************
! rpts: sets up the radial grid
!
! Inputs: ri - dimensional radius of the ICB
!         rs - dimensional radius of the CMB
!
! Outputs: r - array of radius points (global variable)
!
! Ni points are always allocated to the inner core; Nc to outer core
!******************************************************************************
  subroutine rpts_core(ri,rs)

    double precision, intent(in) :: ri, rs
    double precision :: dric, droc               !dr, ic, oc, strat layer
    integer          :: i,bottom,oclower,ocupper !where to start and stop in core

    r_c(0) = 0.d0                                !First point always at the origin

    if(ri .ne. 0.d0) then                        !Check there is an inner core
      dric =     ri /real(Ni)                    !Divide the points equally in space
      do i = 1, Ni
        r_c(i)  = r_c(i-1) + dric
        r2_c(i) = r_c(i)*r_c(i)
        r3_c(i) = r_c(i)*r_c(i)*r_c(i)
        r4_c(i) = r_c(i)*r_c(i)*r_c(i)*r_c(i)
        r5_c(i) = r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)
        r6_c(i) = r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)
        r7_c(i) = r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)
        r8_c(i) = r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)
        r9_c(i) = r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)
        r10_c(i) = r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)
        r11_c(i) = r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)
      enddo
      bottom = Ni+1
    else                                     !If not, OC starts from the centre
      bottom = 1
    endif
                                             !Now do the outer core
    if(ri .ne. 0.d0) then                    !Check there is an inner core
      if(float(ceiling(rs)) == rc) then
        droc = (rs-ri)/real(N  - bottom)
      else
        droc = (rs-ri)/real((Nc+Ni) - bottom)
      endif

      r_c(bottom)  = r_c(Ni)
      r2_c(bottom) = r_c(Ni)*r_c(Ni)
      r3_c(bottom) = r_c(Ni)*r_c(Ni)*r_c(Ni)
      r4_c(bottom) = r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)
      r5_c(bottom) = r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)
      r6_c(bottom) = r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)
      r7_c(bottom) = r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)
      r8_c(bottom) = r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)
      r9_c(bottom) = r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)
      r10_c(bottom) = r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)
      r11_c(bottom) = r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)*r_c(Ni)

      oclower = bottom + 1
    else
      if(float(ceiling(rs)) == rc) then
        droc = (rs-ri)/real(N -bottom+1)
      else
        droc = (rs-ri)/real((Nc+Ni)-bottom+1)
      endif
      oclower = bottom
    endif

    if(float(ceiling(rs)) == rc) then
      ocupper = N
    else
      ocupper = Nc+Ni
    endif

    do i = oclower, ocupper                       !Convecting core
      r_c(i)  = r_c(i-1) + droc
      r2_c(i) = r_c(i)*r_c(i)
      r3_c(i) = r_c(i)*r_c(i)*r_c(i)
      r4_c(i) = r_c(i)*r_c(i)*r_c(i)*r_c(i)
      r5_c(i) = r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)
      r6_c(i) = r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)
      r7_c(i) = r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)
      r8_c(i) = r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)
      r9_c(i) = r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)
      r10_c(i) = r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)
      r11_c(i) = r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)
    enddo
    bottom = ocupper + 1

    if(float(ceiling(rs)) .ne. rc) then
      r_c(bottom)  = r_c(ocupper)
      r2_c(bottom) = r_c(ocupper)*r_c(ocupper)
      r3_c(bottom) = r_c(ocupper)*r_c(ocupper)*r_c(ocupper)
      r4_c(bottom) = r_c(ocupper)*r_c(ocupper)*r_c(ocupper)*r_c(ocupper)
      r5_c(bottom) = r_c(ocupper)*r_c(ocupper)*r_c(ocupper)*r_c(ocupper)*r_c(ocupper)
      r6_c(bottom) = r_c(ocupper)*r_c(ocupper)*r_c(ocupper)*r_c(ocupper)*r_c(ocupper) &
                        *r_c(ocupper)
      r7_c(bottom) = r_c(ocupper)*r_c(ocupper)*r_c(ocupper)*r_c(ocupper)*r_c(ocupper) &
                        *r_c(ocupper)*r_c(ocupper)
      r8_c(bottom) = r_c(ocupper)*r_c(ocupper)*r_c(ocupper)*r_c(ocupper)*r_c(ocupper) &
                        *r_c(ocupper)*r_c(ocupper)*r_c(ocupper)
      r9_c(bottom) = r_c(ocupper)*r_c(ocupper)*r_c(ocupper)*r_c(ocupper)*r_c(ocupper) &
                        *r_c(ocupper)*r_c(ocupper)*r_c(ocupper)*r_c(ocupper)
      r10_c(bottom) = r_c(ocupper)*r_c(ocupper)*r_c(ocupper)*r_c(ocupper)*r_c(ocupper)&
                        *r_c(ocupper)*r_c(ocupper)*r_c(ocupper)*r_c(ocupper)*r_c(ocupper)
      r11_c(bottom) = r_c(ocupper)*r_c(ocupper)*r_c(ocupper)*r_c(ocupper)*r_c(ocupper)&
                        *r_c(ocupper)*r_c(ocupper)*r_c(ocupper)*r_c(ocupper)*r_c(ocupper)*r_c(ocupper)

      droc = (rc-rs)/real(N - bottom)
      do i = bottom+1, N                       !Convecting core
        r_c(i)  = r_c(i-1) + droc
        r2_c(i) = r_c(i)*r_c(i)
        r3_c(i) = r_c(i)*r_c(i)*r_c(i)
        r4_c(i) = r_c(i)*r_c(i)*r_c(i)*r_c(i)
        r5_c(i) = r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)
        r6_c(i) = r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)
        r7_c(i) = r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)
        r8_c(i) = r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)
        r9_c(i) = r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)
        r10_c(i) = r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)
        r11_c(i) = r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)*r_c(i)
      enddo
    endif

  end subroutine rpts_core

!******************************************************************************
! adiabat_poly: Calculates the adiabatic temperature using equation (14) of
!               Davies (PEPI, 2014, submitted)
!
! Inputs: N - index for CMB
!         Tc - Temperature at CMB
!
! Outputs: Qa - Adiabatic heat flux
!          Tar_c - Adiabatic temperature profile (global variable)
!          dTadr_c - Adiabatic temperature gradient (global variable)
!
! The correction procedure for Ta is the same as described in adiabat()
!******************************************************************************
  subroutine adiabat_poly_core(N, Tc, Qa)

    integer         , intent(in) :: N
    double precision, intent(in) :: Tc
    double precision :: Qa
!    double precision :: delta

    T_cen  = Tc / (1.d0 + t1_c*(r_c(N)/1e3) + t2_c*(r2_c(N)/1e6) + t3_c*(r3_c(N)/1e9))
    Tar_c  = T_cen * (1.d0 + t1_c*(r_c/1e3) + t2_c*(r2_c/1e6)    + t3_c*(r3_c/1e9))

!    Tar_c   = T_cen * (1.d0 + t1_c*(r_c/1e3) + t2_c*(r2_c/1e6) + t3_c*(r3_c/1e9))
!    delta = Tc - Tar_c(N)              !Compare Ta @ CMB to Tc

!    T_cen  = T_cen + delta
!    Tar_c   = T_cen * (1.d0 + t1_c*(r_c/1e3) + t2_c*(r2_c/1e6) + t3_c*(r3_c/1e9))

    dTadr_c = T_cen * (t1_c + 2.d0*t2_c*(r_c/1e3) + 3.d0*t3_c*(r2_c/1e6))
    dTadr_c = dTadr_c / 1e3
    Qa    = -4.d0*pi*r2_c(N)*kr_c(N)*dTadr_c(N)

  end subroutine adiabat_poly_core

  subroutine adiabat_ganymede_core(N, Tc, Qa)

    integer         , intent(in) :: N
    double precision, intent(in) :: Tc
    double precision :: Qa
    double precision :: fac(0:N)

    fac     = alphaT_c * gr_c(N) / cp_c
    Tar_c   = Tc * dexp( -alphaT_c * G * rho0_oc * 2d0 * pi * (r2_c - rc**2)/(3d0*cp_c) )
    dTadr_c = Tar_c * fac
    Qa      = -4.d0*pi*r2_c(N)*kr_c(N)*dTadr_c(N)
    write(456,*) Tar_c

  end subroutine adiabat_ganymede_core

  subroutine adiabat_init(Tc, Qa)

    double precision :: Tc, Qa

    T_cen  = Tm_c(Ni) / (1.d0 + t1_c*(r_c(Ni)/1e3) + t2_c*(r2_c(Ni)/1e6) + t3_c*(r3_c(Ni)/1e9))
    Tar_c  = T_cen    * (1.d0 + t1_c*(r_c/1e3)     + t2_c*(r2_c/1e6)     + t3_c*(r3_c/1e9))

    Tc = Tar_c(N)

    dTadr_c = T_cen * (t1_c + 2.d0*t2_c*(r_c/1e3) + 3.d0*t3_c*(r2_c/1e6))
    dTadr_c = dTadr_c / 1e3
    Qa    = -4.d0*pi*r2_c(N)*kr_c(N)*dTadr_c(N)

  end subroutine adiabat_init

!******************************************************************************
! density_poly: Calculate core density using equation (7) of Davies (2014, PEPI)
!
! rhor_c(Ni)*conc*alphac_c is the part of the density jump due to Oxygen.
! Presently this is not used because the density jump input to the code does not
! have to equal that in the density profile.
!******************************************************************************
  subroutine density_poly_core(Ni, conc, alphac)

    integer, intent(in)          :: Ni
    double precision, intent(in) :: conc, alphac
    double precision             :: drho_comp

    rhor_c         = 0.d0

    if(Ni > 0)  then
      rhor_c(0:Ni)   = rho0_ic + rho2_ic*(r2_c(0:Ni)/1e6)
      rhor_c(Ni+1:N) = rho0_oc + rho1_oc*(r_c(Ni+1:N)/1e3) + rho2_oc*(r2_c(Ni+1:N)/1e6) + rho3_oc*(r3_c(Ni+1:N)/1e9)
    else
      rhor_c(0:N)    = rho0_oc + rho1_oc*(r_c(0:N)/1e3)    + rho2_oc*(r2_c(0:N)/1e6)    + rho3_oc*(r3_c(0:N)/1e9)
    endif

    drho_comp = rhor_c(Ni)*conc*alphac

!    write(*,'(A16,E16.6,A16,E16.6)') 'drho_comp=',drho_comp, 'drho_prem=', rhor_c(Ni)-rhor_c(Ni+1)

  end subroutine density_poly_core

!******************************************************************************
! gravity_poly: Calculates gravity using equations (10) and (11) of Davies
!               (2014, PEPI, submitted)
!******************************************************************************
  subroutine gravity_poly_core(Ni)

    integer, intent(in) :: Ni
    double precision    :: M_core(0:N)

    M_core = 0.d0; gr_c=0.d0

    M_core(0:Ni)   = rho0_ic*(r_c(0:Ni)/1e3)/3.d0 + rho2_ic*(r3_c(0:Ni)/1e9)/5.d0

    M_core(Ni+1:N) = M_core(Ni) + &
                    rho0_oc*(r_c (Ni+1:N)/1e3)/3.d0 + rho1_oc*(r2_c(Ni+1:N)/1e6)/4.d0  + &
                    rho2_oc*(r3_c(Ni+1:N)/1e9)/5.d0 + rho3_oc*(r4_c(Ni+1:N)/1e12)/6.d0 - &
                   (rho0_oc*(r_c (Ni+1)  /1e3)/3.d0 + rho1_oc*(r2_c(Ni+1)  /1e6)/4.d0 + &
                    rho2_oc*(r3_c(Ni+1)  /1e9)/5.d0 + rho3_oc*(r4_c(Ni+1)  /1e12)/6.d0)

    gr_c = 4.d0 * pi * G*M_core*1e3

  end subroutine gravity_poly_core

!******************************************************************************
! grav_pot_poly: Calculates the gravitational potential using equation (12) of
!                Davies (2014, PEPI, submitted)
!
! Inputs: Ni - Index for ICB
!         N - Index for CMB
!******************************************************************************
  subroutine grav_pot_poly_core(Ni, N)

    integer, intent(in) :: Ni, N
    double precision    :: gpot_cmb, gpot_icb, gpot_icbi
    integer :: i

    psir_c = 0.d0; gpot_cmb=0.d0; gpot_icb=0.d0; gpot_icbi=0.d0

    gpot_cmb = rho0_oc*(r2_c(N)/1e6)/6.d0    + rho1_oc*(r3_c(N)/1e9)/12.d0 + &
               rho2_oc*(r4_c(N)/1e12)/20.d0  + rho3_oc*(r5_c(N)/1e15)/30.d0

    gpot_icb = rho0_oc*(r2_c(Ni+1)/1e6)/6.d0   + rho1_oc*(r3_c(Ni+1)/1e9)/12.d0 + &
               rho2_oc*(r4_c(Ni+1)/1e12)/20.d0 + rho3_oc*(r5_c(Ni+1)/1e15)/30.d0

    gpot_icbi= rho0_ic*(r2_c(Ni)/1e6)/6.d0     + &
               rho2_ic*(r4_c(Ni)/1e12)/20.d0

    do i = 0, Ni
      psir_c(i) = rho0_ic*(r2_c(i)/1e6)/6.d0 + rho2_ic*(r4_c(i)/1e12)/20.d0 - gpot_cmb + gpot_icb - gpot_icbi
    enddo

    do i = Ni+1, N
      psir_c(i) = rho0_oc*(r2_c(i)/1e6)/6.d0   + rho1_oc*(r3_c(i)/1e9)/12.d0 + &
                rho2_oc*(r4_c(i)/1e12)/20.d0 + rho3_oc*(r5_c(i)/1e15)/30.d0 - gpot_cmb !-gpot_icb + &
!                rho0_ic*(r2_c(Ni)/1e6)/6.d0 + rho2_ic*(r4_c(Ni)/1e12)/20.d0
    enddo

    psir_c = 4.d0 * pi * G * psir_c * 1e6

  end subroutine grav_pot_poly_core

!******************************************************************************
! pressure_poly: Calculates pressure using equation (13) of Davies (2014, PEPI)
!
! Inputs: Ni - Index for ICB
!         N  - Index for CMB
!         Pc - Pressure at centre of Earth
!******************************************************************************
  subroutine pressure_poly_core(Ni, N, Pc)

    double precision, intent(in)  :: Pc
    integer         , intent(in)  :: N, Ni
    double precision              :: one_coeff, two_coeff, three_coeff, four_coeff, five_coeff, six_coeff, seven_coeff
    double precision              :: P_icb, P_cmb, P_icbi

    pr_c = 0.d0

    one_coeff    =       rho0_oc*rho0_oc/6.d0
    two_coeff    = 7.d0 *rho0_oc*rho1_oc/36.d0
    three_coeff  = 2.d0 *rho0_oc*rho2_oc/15.d0  +      rho1_oc*rho1_oc/16.d0
    four_coeff   =       rho0_oc*rho3_oc/10.d0  + 9.d0*rho1_oc*rho2_oc/100.d0
    five_coeff   = 5.d0 *rho1_oc*rho3_oc/72.d0  +      rho2_oc*rho2_oc/30.d0
    six_coeff    = 11.d0*rho2_oc*rho3_oc/210.d0
    seven_coeff  =       rho3_oc*rho3_oc/42.d0

    P_icb        =  4.d6 * pi * G * (one_coeff  *(r2_c(Ni+1)/1e6)  + two_coeff *(r3_c(Ni+1)/1e9)  + &
                   three_coeff*(r4_c(Ni+1)/1e12) + four_coeff *(r5_c(Ni+1)/1e15) + five_coeff*(r6_c(Ni+1)/1e18) + &
                   six_coeff  *(r7_c(Ni+1)/1e21) + seven_coeff*(r8_c(Ni+1)/1e24))

    P_cmb        =  4.d6 * pi * G * (one_coeff  *(r2_c(N)/1e6)      + two_coeff *(r3_c(N)/1e9)     + three_coeff*(r4_c(N)/1e12) + &
                   four_coeff *(r5_c(N)/1e15)    + five_coeff*(r6_c(N)/1e18)    + six_coeff  *(r7_c(N)/1e21)    + &
                   seven_coeff*(r8_c(N)/1e24))

    if(Ni > 0)  then
      pr_c(Ni+1:N)   = Pc + P_cmb -4.d6 * pi * G * &
                     (one_coeff  *(r2_c(Ni+1:N)/1e6)  + two_coeff *(r3_c(Ni+1:N)/1e9)  + three_coeff*(r4_c(Ni+1:N)/1e12) + &
                     four_coeff *(r5_c(Ni+1:N)/1e15) + five_coeff*(r6_c(Ni+1:N)/1e18) + six_coeff  *(r7_c(Ni+1:N)/1e21) + &
                     seven_coeff*(r8_c(Ni+1:N)/1e24))

      one_coeff    =      rho0_ic*rho0_ic/6.d0
      two_coeff    = 8.d0*rho0_ic*rho2_ic/60.d0
      three_coeff  =      rho2_ic*rho2_ic/30.d0

      P_icbi       = 4.d6 * pi * G * (one_coeff*(r2_c(Ni)/1e6) + two_coeff*(r4_c(Ni)/1e12) + three_coeff*(r6_c(Ni)/1e18))

      pr_c(0:Ni)     = Pc + P_cmb - P_icb + P_icbi -4.d6 * pi * G * (one_coeff*(r2_c(0:Ni)/1e6) + two_coeff*(r4_c(0:Ni)/1e12) + &
                                                               three_coeff*(r6_c(0:Ni)/1e18))

    else
      pr_c(0:N)   = Pc + P_cmb -4.d6 * pi * G * &
                     (one_coeff  *(r2_c(0:N)/1e6)  + two_coeff *(r3_c(0:N)/1e9)  + three_coeff*(r4_c(0:N)/1e12) + &
                     four_coeff *(r5_c(0:N)/1e15) + five_coeff*(r6_c(0:N)/1e18) + six_coeff  *(r7_c(0:N)/1e21) + &
                     seven_coeff*(r8_c(0:N)/1e24))
    endif

  end subroutine pressure_poly_core

!******************************************************************************
! mole2massconc: convert the inputted molar concentrations of O, S and Si into
!                mass concentrations that are evolved in time.
!
! Inputs: Molar concentrations (bars; global variables)
! Outputs: Mass concentrations cO, cS, cSi
!
! Mole fraction (#atoms of solute/total #atoms) c' = A'/AL * c, c = mass fraction
!
! Only O is evolved in time because the core chemistry model used here assumes
! that only O is responsible for the compositional part of the ICB density jump
!
! Molecular weights are set as parameters in parameters.f90
!******************************************************************************
  subroutine mole2massconc_core(cbarO_oc, cbarS_oc, cbarSi_oc, cO, cS, cSi)

    double precision, intent(in)  :: cbarO_oc, cbarS_oc, cbarSi_oc
    double precision, intent(out) :: cO, cS, cSi
    double precision              :: Abar

    Abar = cbarO_oc*AO + cbarS_oc*AS + cbarSi_oc*ASi + (1.d0-cbarO_oc-cbarS_oc-cbarSi_oc)*AFe

    cO  = cbarO_oc *((AO*1.d0)/Abar)
    cS  = cbarS_oc *((AS*1.d0)/Abar)
    cSi = cbarSi_oc*((ASi*1.d0)/Abar)

    write(*,'(A16,E16.6)') 'cO=',cO
    write(*,'(A16,E16.6)') 'cS=',cS
    write(*,'(A16,E16.6)') 'cSi=',cSi

  end subroutine mole2massconc_core

  subroutine mass2moleconc_core(cO, cS, cSi, cbarO_oc, cbarS_oc, cbarSi_oc)

    double precision, intent(in)    :: cO, cS, cSi
    double precision, intent(inout) :: cbarO_oc, cbarS_oc, cbarSi_oc
    double precision                :: Abar

    Abar = 1.d0 / (cO/AO + cS/AS + cSi/ASi + (1.d0-cO-cS-cSi)/AFe)

    cbarO_oc  = cO / ((AO*1.d0)/Abar)
    cbarS_oc  = cS / ((AS*1.d0)/Abar)
    cbarSi_oc = cSi /((ASi*1.d0)/Abar)

    write(*,'(A16,E16.6)') 'cbarO=',cbarO_oc
    write(*,'(A16,E16.6)') 'cbarS=',cbarS_oc
    write(*,'(A16,E16.6)') 'cbarSi=',cbarSi_oc

  end subroutine mass2moleconc_core

!******************************************************************************
! mass: Calculates the mass of the inner core (Mic) and outer core (M_oc) using
!       equations (8) and (9) of Davies (2014, PEPI, submitted)
!
! Inputs: Ni - Index for ICB
!         Nc - Index for CMB
!******************************************************************************
  subroutine mass_poly_core(Ni,Ns,N)

    integer, intent(in) :: N,Ns,Ni
    double precision    :: M_oc_t, M_oc_b, M_ic_t, M_s_t, M_s_b

!   Assumes snow and liquid have same density
!   This restricts compositions

    M_s_t  = rho0_oc*(r3_c(N)    /1e9) /3.d0 + rho1_oc*(r4_c(N)   /1e12)/4.d0 + &
             rho2_oc*(r5_c(N)    /1e15)/5.d0 + rho3_oc*(r6_c(N)   /1e18)/6.d0
    M_s_b  = rho0_oc*(r3_c(Ns+1) /1e9) /3.d0 + rho1_oc*(r4_c(Ns+1)/1e12)/4.d0 + &
             rho2_oc*(r5_c(Ns+1) /1e15)/5.d0 + rho3_oc*(r6_c(Ns+1)/1e18)/6.d0
    M_oc_t = rho0_oc*(r3_c(Ns+1) /1e9) /3.d0 + rho1_oc*(r4_c(Ns+1)/1e12)/4.d0 + &
             rho2_oc*(r5_c(Ns+1) /1e15)/5.d0 + rho3_oc*(r6_c(Ns+1)/1e18)/6.d0
    M_oc_b = rho0_oc*(r3_c(Ni+1) /1e9) /3.d0 + rho1_oc*(r4_c(Ni+1)/1e12)/4.d0 + &
             rho2_oc*(r5_c(Ni+1) /1e15)/5.d0 + rho3_oc*(r6_c(Ni+1)/1e18)/6.d0

    M_s    = 4.d0 * pi * (M_s_t  - M_s_b)  * 1e9
    M_oc   = 4.d0 * pi * (M_oc_t - M_oc_b) * 1e9

    M_ic_t = rho0_ic*(r3_c(Ni)/1e9)/3.d0 + rho2_ic*(r5_c(Ni)/1e15)/5.d0
    M_ic   = 4.d0 * pi * (M_ic_t) * 1e9

    M_c    = M_oc + M_ic + M_s

  end subroutine mass_poly_core

!******************************************************************************
! mass_poly_correct: Correct for changes in the core mass over time.
!                 Assume mass of IC is constant.
!
! Inputs: Ni - Index for ICB
!         Nc - Index for CMB
!******************************************************************************

  subroutine mass_poly_correct(ri)

    double precision, intent(in) :: ri
    double precision             :: rho_new, fac, ic, oc1, oc2

    if(ri == ri0) return

    fac =  3.d0 / ( (rc*rc*rc/1e9) - (ri*ri*ri/1e9) )

    oc1 = rho0_oc*( (rc*rc*rc/1e9) - (ri0*ri0*ri0)/1e9)/3.d0

    oc2 = rho1_oc*( (ri0*ri0*ri0*ri0/1e12)         - (ri*ri*ri*ri/1e12) )/4.d0 + &
          rho2_oc*( (ri0*ri0*ri0*ri0*ri0/1e15)     - (ri*ri*ri*ri*ri)/1e15 )/5.d0 + &
          rho3_oc*( (ri0*ri0*ri0*ri0*ri0*ri0/1e18) - (ri*ri*ri*ri*ri*ri)/1e18 )/6.d0

    ic  = rho0_ic*( (ri0*ri0*ri0/1e9)          - (ri*ri*ri)/1e9  )/3.d0 + &
          rho2_ic*( (ri0*ri0*ri0*ri0*ri0/1e15) - (ri*ri*ri*ri*ri)/1e15 )/5.d0

    rho_new = (oc1 - oc2 + ic) * fac
    rho0_oc = rho_new

  end subroutine mass_poly_correct

  subroutine le_mass_test(N, Ns, time, c_liq, c_snow)

    integer         , intent(in)  :: N, Ns
    double precision, intent(in)  :: time, c_snow(0:N), c_liq
    double precision  :: integrand1(0:N), integrand2(0:N)
    integer           :: i

    integrand1=0.d0; integrand2=0.d0
    do i = 0, N
      if(i<Ns+1) integrand1(i) = c_liq * r_c(i)*r_c(i) * rhor_c(i)
      if(i>Ns)   integrand2(i) = c_snow(i) * r_c(i)*r_c(i) * rhor_c(i)
    enddo
    call splat(integrand1,r_c,N+1)
    call splat(integrand2,r_c,N+1)

    write(101,*) time, c_liq, 4.d0 * pi * integrand1(N), 4.d0 * pi * integrand2(N)

  end subroutine le_mass_test

  subroutine solid_mass_test(N, Ns, time, phi, phi_tot)

    integer         , intent(in)  :: N, Ns
    double precision, intent(in)  :: time, phi(0:N)
    double precision  :: integrand1(0:N), integrand2(0:N), phi_tot
    integer           :: i

    integrand1=0.d0; integrand2=0.d0; phi_tot = 0.d0
    do i = 0, N
      if(i>Ns) then 
        integrand1(i) = phi(i) * r_c(i)*r_c(i)
        integrand2(i) =          r_c(i)*r_c(i)
      endif
    enddo
    call splat(integrand1,r_c,N+1)
    call splat(integrand2,r_c,N+1)

    phi_tot = integrand1(N) / integrand2(N)
write(100,*) time, r_c(Ns+1), 4.d0*pi*integrand1(N),  4.d0*pi*integrand2(N), phi_tot

  end subroutine solid_mass_test

  subroutine dMdt(Nsl,Ms,Ml)

    integer, intent(in) :: Nsl
    integer             :: i
    double precision    :: int1(0:N), int2(0:N), int3(0:N), Mc
    double precision, intent(out) :: Ms,Ml

    int1=0.d0; int2=0.d0; int3=0.d0
    do i = 0, N
      if(i.ge.(Nsl+1)) int1(i) = r_c(i)*r_c(i) * rhor_c(i)
      if(i .le. Nsl)   int2(i) = r_c(i)*r_c(i) * rhor_c(i)
      int3(i) = r_c(i)*r_c(i) * rhor_c(i)
    enddo
    call splat(int1,r_c,N)
    Ms=4.d0*pi*int1(N-1)
    call splat(int2,r_c,N)
    Ml=4.d0*pi*int2(N-1)
    call splat(int3,r_c,N)
    Mc=4.d0*pi*int3(N-1)

  end subroutine dMdt

  subroutine dcMdt(Nsl,cliq_s,phi_s,ccore,cMs,cMl)

    double precision, intent(in) :: cliq_s(0:N), phi_s(0:N)
    integer, intent(in) :: Nsl
    double precision, intent(in) :: ccore
    integer             :: i
    double precision    :: int1(0:N), int2(0:N), int3(0:N), ct, cc, cMc
    double precision, intent(out) :: cMs,cMl

    int1=0.d0; int2=0.d0; int3=0.d0
    do i = 0, N
      if(i.ge.(Nsl+1)) then
        cc=0.d0
        ct=cliq_s(i)*(1.d0-phi_s(i))
      else
        cc=ccore
        ct=ccore
      endif

      int1(i) = r_c(i)*r_c(i) * rhor_c(i) * (1.d0-phi_s(i)) * cliq_s(i)
      int2(i) = r_c(i)*r_c(i) * rhor_c(i) * cc ! Liquid core
      int3(i) = r_c(i)*r_c(i) * rhor_c(i) * ct ! Total
    enddo
    call splat(int1,r_c,N)
    cMs=4.d0*pi*int1(N-1)
    call splat(int2,r_c,N)
    cMl=4.d0*pi*int2(N-1)
    call splat(int3,r_c,N)
    cMc=4.d0*pi*int3(N-1)

write(479,*) 'cMs, cMl, cMc = ', cMs, cMl, cMc

  end subroutine dcMdt

!******************************************************************************
! melting: The melting curve using equation (54)
!       of Francis Nimmo's revised treatise article on Energetics of the core
!******************************************************************************

  subroutine melting_core()

    Tm_c    = Tm0_c * (1.d0 + Tm1_c*pr_c + Tm2_c*pr_c*pr_c) - meltdep
    dTmdP_c = Tm0_c * (       Tm1_c   + 2.d0*Tm2_c*pr_c)

  end subroutine melting_core

!******************************************************************************
! melting: The melting curve using DP18
!          DIFFERENT FROM eqn 3 Williams & Nimmo 2004
!******************************************************************************

  subroutine melting_mars(cbar)

    double precision, intent(in) :: cbar
    double precision :: TmFe(0:N)

    Tm_c    = Tm0_c * (1.d0 - ( cbar-meltdep) ) * (1.d0 + Tm1_c*pr_c/1e9 + Tm2_c*pr_c*pr_c/1e18 + Tm3_c*pr_c*pr_c*pr_c/1e27)
    dTmdP_c = Tm0_c * (1.d0 - ( cbar-meltdep) ) * (       Tm1_c     + 2.d0*Tm2_c*pr_c/1e9  + 3.d0*Tm3_c*pr_c*pr_c/1e18)
    dTmdP_c = dTmdP_c / 1e9

  end subroutine melting_mars

  subroutine melting_ganymede(cbar)

    double precision, intent(in) :: cbar
    double precision :: A(0:N), B(0:N), C(0:N), D(0:N), E(0:N)
    double precision :: dA(0:N), dB(0:N), dC(0:N), dD(0:N), dE(0:N)

    A = -2.4724*(pr_c/1e9)**4 + 28.025*(pr_c/1e9)**3 + 9.1404*(pr_c/1e9)**2 + 581.71*(pr_c/1e9) + 3394.8
    B =  1.7978*(pr_c/1e9)**4 - 6.7881*(pr_c/1e9)**3 - 197.69*(pr_c/1e9)**2 - 271.69*(pr_c/1e9) - 8219.5
    C = -0.1702*(pr_c/1e9)**4 - 9.3959*(pr_c/1e9)**3 + 163.53*(pr_c/1e9)**2 - 319.35*(pr_c/1e9) + 5698.6
    D = -0.2308*(pr_c/1e9)**4 + 7.1000*(pr_c/1e9)**3 - 64.118*(pr_c/1e9)**2 + 105.98*(pr_c/1e9) - 1621.9
    E =  0.2302*(pr_c/1e9)**4 - 5.3688*(pr_c/1e9)**3 + 38.124*(pr_c/1e9)**2 - 46.681*(pr_c/1e9) + 1813.8
    Tm_c    = A*cbar**4 + B*cbar**3 + C*cbar**2 + D*cbar + E

    dA = -4.0*2.4724*(pr_c/1e9)**3 + 3.0*28.025*(pr_c/1e9)**2 + 2.0*9.1404*(pr_c/1e9) + 581.71
    dB =  4.0*1.7978*(pr_c/1e9)**3 - 3.0*6.7881*(pr_c/1e9)**2 - 2.0*197.69*(pr_c/1e9) - 271.69
    dC = -4.0*0.1702*(pr_c/1e9)**3 - 3.0*9.3959*(pr_c/1e9)**2 + 2.0*163.53*(pr_c/1e9) - 319.35
    dD = -4.0*0.2308*(pr_c/1e9)**3 + 3.0*7.1000*(pr_c/1e9)**2 - 2.0*64.118*(pr_c/1e9) + 105.98
    dE =  4.0*0.2302*(pr_c/1e9)**3 - 3.0*5.3688*(pr_c/1e9)**2 + 2.0*38.124*(pr_c/1e9) - 46.681

    dTmdP_c = dA*cbar**4 + dB*cbar**3 + dC*cbar**2 + dD*cbar + dE

    dTmdP_c = dTmdP_c / 1e9

  end subroutine melting_ganymede

  subroutine melting_mars_snowzone(Nsnow)

    integer, intent(in) :: Nsnow

    Tm_c(Nsnow+1:N) = Tar_c(Nsnow+1:N)

  end subroutine melting_mars_snowzone

!******************************************************************************
! melting_poly: The melting curve using equation (16) of Davies (2014, submitted
!               to PEPI)
!******************************************************************************
  subroutine entropy_melting()

    ds_c = 0.d0

    ds_c = ds0_c + ds1_c*(pr_c/1e9) + ds2_c*(pr_c*pr_c)/1e18 + ds3_c*(pr_c*pr_c*pr_c)/1e27 + ds4_c*(pr_c*pr_c*pr_c*pr_c)/1e36

  end subroutine entropy_melting

!******************************************************************************
! melting_poly: The melting curve using equation (16) of Davies (2014, submitted
!               to PEPI)
!******************************************************************************
  subroutine melting_poly_core()

    Tm_c = 0.d0; dTmdP_c = 0.d0

    Tm_c    = Tm0_c + Tm1_c*(pr_c/1e9) + Tm2_c*(pr_c*pr_c)/1e18 + Tm3_c*(pr_c*pr_c*pr_c)/1e27 + Tm4_c*(pr_c*pr_c*pr_c*pr_c)/1e36
    dTmdP_c =         Tm1_c            + 2.d0*Tm2_c*pr_c/1e9    + 3.d0*Tm3_c*pr_c*pr_c/1e18   + 4.d0*Tm4_c*(pr_c*pr_c*pr_c)/1e27
    dTmdP_c = dTmdP_c / 1e9
    Tm_c    = Tm_c - meltdep

  end subroutine melting_poly_core

!******************************************************************************
! Solid_conc: get cS from cL
!******************************************************************************
  subroutine solid_conc(cl, lambdal, lambdas, dmu, cs)

    double precision, intent(in)  :: cl, lambdal, lambdas, dmu
    double precision, intent(out) :: cs
    double precision              :: xmin, xmax, TmFe, dsri

    xmin = 0.1d-10
    xmax = 1.2d0
    TmFe = Tm_c(Ni+1)
    dsri = ds_c(Ni+1)

    if(ri0 < 1e-10) then
      TmFe = Tm_c(0)
      dsri = ds_c(0)
    endif

    cs = 0.d0
    if(cl==0.d0) then
      cs = 0.d0
      return
    endif

    cs = rtbis(xmin, xmax, 1d-10, lambdal, lambdas, dmu ,cl ,TmFe, dsri)

  end subroutine solid_conc

  double precision function f(x,lambdal,lambdas,dmu,cl,TmFe,dSFe)
    implicit none

    double precision, intent(in) :: x, TmFe, dSFe, cl, lambdas, lambdal, dmu

    f  = dmu + lambdal*cl - lambdas*x  - kb*TmFe*dlog(x/cl)*( dSFe + (x-cl) ) /dSFe

  end function f

  double precision function rtbis( x1, x2, xacc, lambdal,lambdas,dmu,cl,TmFe,dSFe)

    integer , parameter :: n=200
    integer             :: i
    double precision    :: x1, x2, xacc
    double precision    :: fst, fmid, dx, xmid
    double precision, intent(in) :: TmFe, dSFe, cl, lambdas, lambdal, dmu

    fmid = f(x2,lambdal,lambdas,dmu,cl,TmFe,dSFe)
    fst  = f(x1,lambdal,lambdas,dmu,cl,TmFe,dSFe)

    if(fst*fmid .ge. 0.d0) stop 'Root must be bracketed for bisection'

    if(fst .le. 0.d0) then
      rtbis = x1
      dx   = x2-x1
    else
      rtbis = x2
     dx   = x1-x2
    endif

    do i = 1, n
      dx = dx * 0.5d0
      xmid = rtbis + dx
      fmid = f(xmid,lambdal,lambdas,dmu,cl,TmFe,dSFe)

      if(fmid .le. 0.d0) rtbis=xmid
      if(dabs(dx) .lt. xacc .or. fmid .eq. 0.d0) return
    enddo

  end function rtbis

!******************************************************************************
! melting_pt_dep: The melting point depression for each light element
!******************************************************************************
  subroutine melting_pt_dep(cl, cs, dTm)

    double precision, intent(in)  :: cs, cl
    double precision, intent(out) :: dTm(0:N)

    dTm = 0.d0
    dTm = (Tm_c/ds_c) *  (cs - cl)

  end subroutine melting_pt_dep

!******************************************************************************
! ricalc: find the position  of the ICB from the intersection of Ta and Tm
!
! Inputs: N - total number of gridpoints
!
! Outputs: ri - ICB radius
!******************************************************************************
  subroutine ricalc(N,ri)

    implicit none

    integer, intent(in)           :: N
    double precision, intent(out) :: ri
    double precision              :: Tdiff
    integer                       :: i, ript

    Tdiff  = 0.d0

!   The two T profiles are on the same grid at each time pt
!   so we can just compare then pointwise

    ript = 0
    a: do i = 0, N
         Tdiff = Tar_c(i) - Tm_c(i)
         if(Tdiff .ge. 0.d0 ) then
           ript   = i
           exit a
         endif
       enddo a

    ri = r_c(ript)

!   ceiling rounds up. Was getting ri = rc-1e-7!
    if(ceiling(ri) .ge. rc) stop "ri > rc!"

    write(*,'(A16,E16.6,I16)') 'ri =',ri, ript

  end subroutine ricalc

! Find root of Ta-Tm = 0
! Approximate Tm as locally linear
  double precision function Td(x, xm0, mm, cm, xa0, ma, ca)
    implicit none

    double precision, intent(in) :: x, xm0, mm, cm, xa0, ma, ca

    Td = ( (ma*x)-(ma*xa0) + ca )  - ( (mm*x)-(mm*xm0) + cm )

  end function Td

  double precision function Tdbis( x1, x2, xacc, xm0, mm, cm, xa0, ma, ca)

    integer , parameter :: n=200
    integer             :: i
    double precision    :: x1, x2, xacc, x, xm0, mm, cm, xa0, ma, ca
    double precision    :: fst, fmid, dx, xmid

    fmid = Td(x2, xm0, mm, cm, xa0, ma, ca)
    fst  = Td(x1, xm0, mm, cm, xa0, ma, ca)

    if(fst*fmid .ge. 0.d0) then 
      print *, 'Root not bracketed', fst, fmid
      return
    endif

    if(fst .le. 0.d0) then
      Tdbis = x1
      dx   = x2-x1
    else
      Tdbis = x2
      dx   = x1-x2
    endif

    do i = 1, n
      dx = dx * 0.5d0
      xmid = Tdbis + dx
      fmid = Td(xmid, xm0, mm, cm, xa0, ma, ca)

      if(fmid .le. 0.d0) Tdbis=xmid
      if(dabs(dx) .lt. xacc .or. fmid .eq. 0.d0) return
    enddo

  end function Tdbis


!******************************************************************************
! rscalc: find the position  of the ICB from the intersection of Ta and Tm
!
! Inputs: N - total number of gridpoints
!
! Outputs: rs - Outer Boundary radius
!******************************************************************************
  subroutine rscalc(N, rt, rm1, rm2, rb)

    implicit none

    integer, intent(in)           :: N
    double precision, intent(out) :: rt, rm1, rm2, rb
    double precision              :: Tdiff, Tmgrad, Tagrad, mT
    integer                       :: rspt1, rspt2, rspt3, rspt4, i
    integer                       :: rspta, rsptb, rsptc, rsptd

!   Check that we are not liquid at the top
!   In 'Earth' case this would mean whole core is solid
!    if( Tm_c(N) > Tar_c(N) ) return

    Tdiff  = 0.d0; mT=-1.d0

!   The two T profiles are on the same grid at each time pt
!   so we can just compare then pointwise
    rspt1=0; rspt2=0; rspt3=0; rspt4=0
    rspta=0; rsptb=0; rsptc=0; rsptd=0
    do i = N, 0, -1                              !Start at CMB and work down
      Tdiff = Tar_c(i) - Tm_c(i)
      if(Tdiff > mT) mT = Tdiff
      if( (Tdiff > 0.d0) .and. (rspta==0) ) then !Intersect is between rspt1 and rspt1+1
        rspt1   = i
        rspta   = 1
      endif
      if( (Tdiff < 0.d0) .and. (rspta==1) .and. (rsptb==0) ) then
        rspt2   = i
        rsptb   = 1
      endif
    enddo

    do i = 0, N                                   !Start at centre
      Tdiff = Tar_c(i) - Tm_c(i)
      if( (Tdiff > 0.d0) .and. (rsptc==0) ) then  !If liquid and not found liquid before
        rspt3   = i                               !Solid inner core scenario
        rsptc   = 1                               !...and stop looking
      endif
      if( (Tdiff < 0.d0) .and. (rsptc==1) .and. (rsptd==0) ) then !If solid, not found solid before, and previously found liquid
        rspt4   = i                               !Base of middle solid layer
        rsptd   = 1                               !...and stop looking
      endif
    enddo

    write(*,'(A16,2E16.6)') 'rs (naive search) = ', r_c(rspt1), rc
    Tagrad = ( Tar_c(rspt1+1)-Tar_c(rspt1) ) / ( r_c(rspt1+1) - r_c(rspt1) )
    Tmgrad = ( Tm_c(rspt1+1) -Tm_c(rspt1)  ) / ( r_c(rspt1+1) - r_c(rspt1) )
    rt = Tdbis( r_c(rspt1), r_c(rspt1+1), 1d-5, r_c(rspt1), Tmgrad, Tm_c(rspt1), r_c(rspt1), Tagrad, Tar_c(rspt1) )

!   Protects against predicting a layer when its just a rounding problem. 
    if(float(ceiling(r_c(rspt1))) == rc) rt=rc

    rm1 = r_c(rspt2)
    rm2 = r_c(rspt4-1)
    rb  = r_c(rspt3)

    if( (ri0 > 0.d0) .and. (mT < 0.d0) ) rb = rc !Whole core is below the melting pt
                                                 !Using ri0 shows there was an IC at a prev time

  end subroutine rscalc

!******************************************************************************
! crfac: Calculate Cr using equation (50)
!       of Francis Nimmo's revised treatise article on Energetics of the core
!
! Inputs: Ni - Index for ICB radius
!         Tc - CMb temperature
!         ri - Current ICB radius
!
! Outputs: Cr
!******************************************************************************
  subroutine crfac(Ni,Tc,ri,Cr)

    integer, intent(in)           :: Ni
    double precision, intent(in)  :: Tc,ri
    double precision, intent(out) :: Cr
    double precision :: dTadP

! 
!   OK to do integral from 0 to N as dcldT_s=0 outside snow zone  If the temperature at the centre of the Earth is not below the
!   melting temperature at the same pressure then we have no IC
    if(ri==0.d0 .or. (float(ceiling(ri))==rc)) then
      dTadP = 0.d0
      Cr = 0.d0
    else
      dTadP = -dTadr_c(Ni+1)/(rhor_c(Ni+1)*gr_c(Ni+1)) !Adiabat at ICB
      Cr    = -(Tar_c(Ni+1)/Tc)*(1.d0/(rhor_c(Ni+1)*gr_c(Ni+1)))*(1.d0/(dTmdP_c(Ni+1) - dTadP))
    endif

    write(*,'(A16,E16.6)') 'Cr =',Cr

    write(14,'(11E16.6)') r_c(Ni), dTadP, dTmdP_c(Ni+1),rhor_c(Ni+1),Tar_c(Ni+1),gr_c(Ni+1),pr_c(Ni+1),Cr,dTadr_c(Ni+1), &
                         -dTmdP_c(Ni+1)*rhor_c(Ni+1)*gr_c(Ni+1), T_cen

  end subroutine crfac

!******************************************************************************
! ccfac: Calculate Cc using equation (51)
!        of Francis Nimmo's revised treatise article on Energetics of the core
!
! Inputs: Ni - Index for ICB
!         conc - O concentration
!
! Outputs: Cc
!******************************************************************************
  subroutine ccfac(Ni, cl, cs, Cc)

    integer, intent(in) :: Ni
    double precision, intent(in)  :: cl, cs
    double precision, intent(out) :: Cc

    Cc = 0.d0
    Cc = 4.d0*pi*r2_c(Ni+1)*rhor_c(Ni+1) * (cl  - cs) / M_oc

    write(*,'(A16,E16.6)') 'Cc =',Cc

  end subroutine ccfac

!******************************************************************************
! cpfac: Calculate integral of rho*xi/(xl**2) d xl / dt * (Tcmb/Ta)
!
! Inputs: Ni - Index for ICB
!         conc - O concentration
!
! Outputs: Cc
!******************************************************************************
  subroutine cpfac(Nsl, cliq, dcldT_s, Tc, Cp)

    integer, intent(in) :: Nsl
    integer             :: i
    double precision    :: int1(0:N), int2(0:N)
    double precision, intent(in)  :: cliq(0:N), Tc, dcldT_s(0:N)
    double precision, intent(out) :: Cp

!   Chk'd int rho r^2 on 5th Dec 15
!   OK to do integral from 0 to N as dcldT_s=0 outside snow zone
    Cp = 0.d0
    int1=0.d0
    do i = 0, N
      int1(i) = ( r_c(i)*r_c(i)*rhor_c(i)*dcldT_s(i) * Tar_c(i) )
    enddo
    call splat(int1,r_c,N+1)

    Cp = -( 4.d0 * pi * int1(N) ) / ( M_oc * Tc )

  end subroutine cpfac

  subroutine liquidusfac(N,Ns,L,xl,phi,fac)

    integer         , intent(in)  :: N, Ns
    double precision, intent(in)  :: L(0:N), xl(0:N), phi(0:N)
    double precision, intent(out) :: fac(0:N)
    double precision              :: dmubar_dcbarS(0:N), dmudxi(0:N)

    dmubar_dcbarS = kb*Tar_c/xl
    dmudxi        = dmubar_dcbarS  * Ev * Na * (1000.d0/AS)

    fac(0:Ns)   = 0.d0
    fac(Ns+1:N) = -(L(Ns+1:N) * (1.d0-phi(Ns+1:N)) ) / ( xl(Ns+1:N) * dmudxi(Ns+1:N) * Tar_c(Ns+1:N) )

  end subroutine liquidusfac

  subroutine xliq_snow(Nsl, Tar_c, cliq_s, dcldT_s)

    integer, intent(in)           :: Nsl
    double precision, intent(in)  :: Tar_c(0:N)
    double precision, intent(out) :: cliq_s(0:N), dcldT_s(0:N)
    double precision              ::  minc, maxc, Tm_Fe(0:N)
    integer :: i

    Tm_Fe = Tm0_c * (1.d0 + Tm1_c*pr_c/1e9 + Tm2_c*pr_c*pr_c/1e18 + Tm3_c*pr_c*pr_c*pr_c/1e27)

    cliq_s(0:Nsl)   = 0.d0
    cliq_s(Nsl+1:N) = 1.d0 + meltdep - ( Tar_c(Nsl:N)/Tm_Fe(Nsl:N) )

    minc = minval(cliq_s)
    maxc = maxval(cliq_s)

    if( (minc < 0.0) .or. (maxc > 0.18) ) then
      write(*,*)  minval(Tar_c(Nsl+1:N)), maxval(Tar_c(Nsl+1:N))
      write(*,*) 'WARNING: Min cl = ', minc, ' ; Max cl = ', maxc
!      stop "Conc out of range"
    endif

    dcldT_s(0:Nsl)   = 0.d0
    dcldT_s(Nsl+1:N) = -1.d0 / ( Tm0_c*(1.d0 + Tm1_c*pr_c(Nsl+1:N)/1e9 + Tm2_c*pr_c(Nsl+1:N)*pr_c(Nsl+1:N)/1e18) + & 
                                               Tm3_c*pr_c(Nsl+1:N)*pr_c(Nsl+1:N)*pr_c(Nsl+1:N)/1e27 ) 

  end subroutine xliq_snow

! x = (1 - phi)*xl
! Since all solid is assumed to fall out of the the snow layer after
! each timestep we take x0 to be the conc from the previous timestep
  subroutine phi_snow(Nsl, Tc, dTcdt, cS, cliq_s, c0_s, phi_s)

    double precision, intent(in) :: Tc, dTcdt, cliq_s(0:N), cS
    integer         , intent(in) :: Nsl
    double precision, intent(out) :: c0_s(0:N), phi_s(0:N)
    double precision :: T_cen, Tar_c_old(0:N), Abar(0:N), dc
    integer :: i

!   Get T from previous timestep
    if( (float(ceiling(rs)) .ne. rc) ) then
      Tar_c_old = 0.d0
      T_cen     = (Tc- dTcdt*secinyr*dt) / (1.d0 + t1_c*(r_c(N)/1e3) + t2_c*(r2_c(N)/1e6) + t3_c*(r3_c(N)/1e9))
      Tar_c_old = T_cen * (1.d0 + t1_c*(r_c/1e3) + t2_c*(r2_c/1e6)    + t3_c*(r3_c/1e9))
      call xliq_snow(Nsl, Tar_c_old, c0_s, Tar_c_old)
    else
      c0_s = cS
    endif

    phi_s(0:Nsl)   = 0.d0
    phi_s(Nsl+1:N) = 1.d0 - (c0_s(Nsl+1:N) / cliq_s(Nsl+1:N))

  end subroutine phi_snow

!******************************************************************************
! dcdt: Calculate the time change of concentration, dcdt, using equation (51)
!       of Francis Nimmo's revised treatise article on Energetics of the core
!
! Inputs: Cc - (see ccfac)
!         Cr - (see crfac)
!         dTcdt - CMB cooling rate
!         conc - O ocncentration
!         forback - Going forwards (=1) or backwards (=2) in time
!******************************************************************************
  subroutine dcdt(Cc, Cr, Cp, dTcdt, conc, forback)

    double precision, intent(in)  :: Cc, Cr, Cp, dTcdt
    double precision, intent(out) :: conc
    double precision              :: dconcdt
    integer                       :: forback

    dconcdt = (Cp + Cc*Cr) * dTcdt * secinyr

    if(forback==1) conc    = conc + dconcdt * dt
    if(forback==2) conc    = conc - dconcdt * dt

    write(*,'(A16, E16.6)') 'Dc/Dt = ', dconcdt
  end subroutine dcdt

!******************************************************************************
! drhodt: Calculate time change of density. Assumes thermal and pressure
!         contributions are negligible
!******************************************************************************
  subroutine drhodt(dconcdt,alphac)

    double precision, intent(in) :: dconcdt, alphac
    double precision             :: ddensitydt

    ddensitydt = dconcdt * rhor_c(N) * alphac

    write(14,'(A16, E16.6)') 'D(rho)/Dt = ', ddensitydt

  end subroutine drhodt

  subroutine secular_num(N,Tc,Qs,Es)

    integer, intent(in) :: N
    double precision, intent(in)  :: Tc
    double precision, intent(out) :: Qs, Es
    double precision  :: integrand1(0:N), integrand2(0:N)
    integer           :: i

    Qs=0.d0; Es=0.d0

    integrand1=0.d0; integrand2=0.d0
    do i = 0, N
      integrand1(i) = cp_c * r_c(i)*r_c(i) * rhor_c(i)*Tar_c(i)
      integrand2(i) = cp_c * r_c(i)*r_c(i) * rhor_c(i) * (Tar_c(i)/Tc - 1.d0)
    enddo
    call splat(integrand1,r_c,N+1)
    call splat(integrand2,r_c,N+1)

    Qs = -4.d0 * pi * integrand1(N) / Tc
    Es = -4.d0 * pi * integrand2(N) / Tc !cp_c * (M_c - rhoT/(Tc*cp_c)) /Tc

  end subroutine secular_num

  subroutine secular_snow(N,Ns,Tc,L,xl,lfac,Qs,Es)

    integer, intent(in)           :: N, Ns
    double precision, intent(in)  :: Tc, L(0:N), lfac(0:N), xl(0:N)
    double precision, intent(out) :: Qs, Es
    double precision  :: integrand1(0:N), integrand2(0:N), integrand3(0:N)
    double precision  :: cp_snow(0:N)
    integer           :: i, Abar

    Qs=0.d0; Es=0.d0

    cp_snow(0:Ns)   = 0.d0
    cp_snow(Ns+1:N) = -L(Ns+1:N)*lfac(Ns+1:N) / xl(Ns+1:N)

    integrand1=0.d0; integrand2=0.d0; integrand3=0.d0
    do i = 0, N
      integrand1(i) = cp_snow(i) * r_c(i)*r_c(i) * rhor_c(i) *  Tar_c(i)
      integrand2(i) = cp_snow(i) * r_c(i)*r_c(i) * rhor_c(i) * (Tar_c(i)/Tc - 1.d0)
    enddo
    call splat(integrand1,r_c,N+1)
    call splat(integrand2,r_c,N+1)

    Qs = -4.d0 * pi * ( integrand1(N) ) / Tc
    Es = -4.d0 * pi * ( integrand2(N) ) / Tc !cp_c * (M_c - rhoT/(Tc*cp_c)) /Tc

  end subroutine secular_snow

  subroutine secular_poly_core(N,Ni,Tc,Qs,Es)

    integer, intent(in)           :: N, Ni
    double precision, intent(in)  :: Tc
    double precision, intent(out) :: Qs, Es
    double precision              :: one_coeff, two_coeff, three_coeff, four_coeff, five_coeff, six_coeff, seven_coeff
    double precision              :: Soc_cmb, Soc_icb, Sic_icb, int_rhoT_dV

    Qs=0.d0; Es=0.d0

    one_coeff=0.d0; two_coeff=0.d0; three_coeff=0.d0; four_coeff=0.d0; five_coeff=0.d0; six_coeff=0.d0; seven_coeff=0.d0

    one_coeff    = rho0_oc*T_cen
    two_coeff    = rho0_oc*T_cen*t1_c   + rho1_oc*T_cen
    three_coeff  = rho2_oc*T_cen      + rho1_oc*T_cen*t1_c   + rho0_oc*T_cen*t2_c
    four_coeff   = rho2_oc*T_cen*t1_c   + rho1_oc*T_cen*t2_c   + rho3_oc*T_cen    + rho0_oc*T_cen*t3_c
    five_coeff   = rho2_oc*T_cen*t2_c   + rho3_oc*T_cen*t1_c   + rho1_oc*T_cen*t3_c
    six_coeff    = rho3_oc*T_cen*t2_c   + rho2_oc*T_cen*t3_c
    seven_coeff  = rho3_oc*T_cen*t3_c

    Soc_cmb = one_coeff*(r3_c(N)/1e9)/3.d0  +  two_coeff*(r4_c(N)/1e12)/4.d0 + three_coeff*(r5_c(N)/1e15)/5.d0 + &
              four_coeff*(r6_c(N)/1e18)/6.d0 + five_coeff*(r7_c(N)/1e21)/7.d0 + six_coeff*(r8_c(N)/1e24)/8.d0 + &
              seven_coeff*(r9_c(N)/1e27)/9.d0
    Soc_icb = one_coeff*(r3_c(Ni)/1e9)/3.d0  +  two_coeff*(r4_c(Ni)/1e12)/4.d0 + three_coeff*(r5_c(Ni)/1e15)/5.d0 + &
              four_coeff*(r6_c(Ni)/1e18)/6.d0 + five_coeff*(r7_c(Ni)/1e21)/7.d0 + six_coeff*(r8_c(Ni)/1e24)/8.d0 + &
              seven_coeff*(r9_c(Ni)/1e27)/9.d0

    one_coeff=0.d0; two_coeff=0.d0; three_coeff=0.d0; four_coeff=0.d0; five_coeff=0.d0; six_coeff=0.d0; seven_coeff=0.d0

    one_coeff    = rho0_ic*T_cen
    two_coeff    = rho0_ic*T_cen*t1_c
    three_coeff  = rho2_ic*T_cen      + rho0_ic*T_cen*t2_c
    four_coeff   = rho2_ic*T_cen*t1_c   + rho0_ic*T_cen*t3_c
    five_coeff   = rho2_ic*T_cen*t2_c
    six_coeff    = rho2_ic*T_cen*t3_c

    Sic_icb = one_coeff*(r3_c(Ni)/1e9)/3.d0  +  two_coeff*(r4_c(Ni)/1e12)/4.d0 + three_coeff*(r5_c(Ni)/1e15)/5.d0 + &
              four_coeff*(r6_c(Ni)/1e18)/6.d0 + five_coeff*(r7_c(Ni)/1e21)/7.d0 + six_coeff*(r8_c(Ni)/1e24)/8.d0

    int_rhoT_dV = 4.d0 * pi * (Soc_cmb - Soc_icb + Sic_icb) * 1e9

    Qs = -cp_c * int_rhoT_dV / Tc

    Es =  cp_c * (M_c - int_rhoT_dV/Tc) /Tc

  end subroutine secular_poly_core

  subroutine radiogenic_decay_core(time, h0_c, h)

    double precision, intent(in)    :: time, h0_c
    double precision, intent(inout) :: h
    double precision                :: t

    t = (ttot/1e6) - time

    h = h0_c * 2**(t/hl)

  end subroutine radiogenic_decay_core

  subroutine radiogenic_decay2(time, h0, h)

    double precision, intent(in)    :: time, h0
    double precision, intent(inout) :: h

    h = h0 * 2**(-time/hl)

  end subroutine radiogenic_decay2

  subroutine radiogenic_num_core(N,Tc,h,Qr,Er)

    integer, intent(in)           :: N
    double precision, intent(in)  :: Tc, h
    double precision, intent(out) :: Qr, Er
    double precision              :: integrand(0:N)
    integer :: i

    Qr = M_c*h

    do i = 0, N
      integrand(i) = r_c(i)*r_c(i)*rhor_c(i)*(1.d0/Tar_c(N) - 1.d0/Tar_c(i))
    enddo
    call splat(integrand,r_c,N)
    Er=4.d0*pi*integrand(N-1)*h!1.d12

  end subroutine radiogenic_num_core

  subroutine lh_coeff_core(N,L)

    integer         , intent(in)  :: N
    double precision, intent(out) :: L(0:N)
    integer :: i

    L = 0.d0
    L = Tm_c * ds_c * kB * Ev * Na * 1000.d0  / AFe

  end subroutine lh_coeff_core

  subroutine latent_core(Ni,L,Tc,dridt,Ql,El)

    integer         , intent(in)  :: Ni
    double precision, intent(in)  :: Tc, dridt, L(0:N)
    double precision, intent(out) :: Ql, El

    Ql = 4.d0*pi*r2_c(Ni+1)*L(Ni+1)*rhor_c(Ni+1)*dridt

    El = Ql*(Tar_c(Ni+1)-Tc)/(Tar_c(Ni+1)*Tc)

  end subroutine latent_core

  subroutine gravitational_freezing(N,Ni,Tc,Cr,cliq,dcldT_s,alphac,Qg,Eg)

    integer, intent(in)          :: N, Ni
    double precision, intent(in) :: Tc, alphac, Cr
    double precision             :: integrand(0:N), dcldT_s(0:N), cliq(0:N)
    double precision             :: Qg_b, Qg_i
    double precision, intent(out):: Qg, Eg
    integer                      :: i
 
    integrand=0.d0
    do i = 0, N
      integrand(i)=psir_c(i)*rhor_c(i)*dcldT_s(i)*Tar_c(i)*r_c(i)*r_c(i)
    enddo
    call splat(integrand,r_c,N+1)
    Qg_i = 4.d0 * pi * (integrand(N)-integrand(Ni+1)) / Tc
!    Qg_b = 4.d0 * pi * r_c(Ni+1)*r_c(Ni+1) * rhor_c(Ni+1) * psir_c(Ni+1) * cliq(Ni+1) * Cr

    ! CHK - Signs here. 
    ! I think its ok...for the liquid region the integral is over the volume MINUS 
    ! the shell that froze into the liquid
    ! for the snow zone the integral is over the snow zone PLUS the shell that 
    ! froze, expanding the snow zone. 
    Qg = -(Qg_i) * alphac
!    Qg = (Qg_i + Qg_b) * alphac
    Eg = Qg/Tc

    write(*,'(A16,3E16.6)') 'Qgs freeze = ', Qg, Qg_i*alphac, Qg_b*alphac

    !**************************************************************
    ! I ran this test to see how to do a volume integral with splat
    ! The answer is to do integrand2(N)-integrand2(Ni+1)
    integrand=0.d0
    do i = 0, N
      integrand(i) = r_c(i) * r_c(i)
    enddo
    call splat(integrand,r_c,N+1)
    write(*,'(A16,3E16.6)') 'CHK vol int: ', integrand(N)-integrand(Ni+1), (r_c(n)**3-r_c(Ni+1)**3)/3.d0
    !**************************************************************

  end subroutine gravitational_freezing

  subroutine gravitational_remelting(N,Ni,Tc,Cr,Cc,Cp,c,alphac,Qg,Eg)

    integer, intent(in)          :: N, Ni
    double precision, intent(in) :: Tc, Cr, Cc, Cp, alphac, c
    double precision             :: integrand(0:N), Qg_i, Qg_b
    double precision, intent(out):: Qg, Eg
    integer                      :: i

    do i = 0, N
      integrand(i)=psir_c(i)*rhor_c(i)*r_c(i)*r_c(i)
    enddo
    call splat(integrand,r_c,N+1)

    ! Use the first form of Get04 eqn 38 as eqn 9 does not hold for slurry
    Qg_i = 4.d0 * pi * (integrand(N)-integrand(Ni+1)) * (Cr*Cc + Cp)
    Qg_b = 4.d0 * pi * r_c(N)*r_c(N) * rhor_c(N) * psir_c(N) * c * Cr

    Qg = (Qg_i - Qg_b) * alphac
    Eg = Qg/Tc

    write(*,'(A16,3E16.6)') 'Qgs remelt = ', Qg, Qg_i*alphac, -Qg_b*alphac

  end subroutine gravitational_remelting

  subroutine gravitational_liq(N,Ni,Tc,Cr,Cc,Cp,alphac,Qg,Eg)

    integer, intent(in)          :: N, Ni
    double precision, intent(in) :: Tc, Cr, Cc, Cp, alphac
    double precision             :: integrand(0:N)
    double precision, intent(out):: Qg, Eg
    integer                      :: i

    do i = 0, N
      integrand(i)=psir_c(i)*rhor_c(i)*r_c(i)*r_c(i)
    enddo
    call splat(integrand,r_c,N+1)
    Qg = 4.d0 * pi * (integrand(N)-integrand(Ni+1)) 
    Qg = (Qg - M_oc * psir_c(Ni)) * (Cr*Cc + Cp) * alphac
    Eg = Qg/Tc

  end subroutine gravitational_liq

! ASSUMES - Only O contributes to the density jump.
  subroutine gravitational_poly_core(N, Ni, Tc, Cr, Cc, Cp, alphac, Qg, Eg)

    integer, intent(in)           :: N,Ni
    double precision, intent(in)  :: Tc, Cr, Cc, alphac, Cp
    double precision, intent(out) :: Qg, Eg

    double precision    :: gpot_cmb,rho_psicmb_cmb,rho_psicmb_icb,rho_psi_cmb, rho_psi_icb
    double precision    :: five_coeff,six_coeff,seven_coeff,eight_coeff,nine_coeff,ten_coeff,eleven_coeff
    double precision    :: int_rhopsi_dV

    gpot_cmb=0.d0; rho_psicmb_cmb=0.d0; rho_psicmb_icb=0.d0; rho_psi_cmb=0.d0; rho_psi_icb=0.d0

!    gpot_cmb = -rho0_oc*(r2_c(N)/1e6)/6.d0   - rho1_oc*(r3_c(N)/1e9)/12.d0 - &
!                rho2_oc*(r4_c(N)/1e12)/20.d0 - rho3_oc*(r5_c(N)/1e15)/30.d0

!    gpot_icb = rho0_oc*(r2_c(Ni+1)/1e6)/6.d0   + rho1_oc*(r3_c(Ni+1)/1e9)/12.d0 + &
!               rho2_oc*(r4_c(Ni+1)/1e12)/20.d0 + rho3_oc*(r5_c(Ni+1)/1e15)/30.d0

!    gpot_icbi= rho0_ic*(r2_c(Ni  )/1e6)/6.d0   + &
!               rho2_ic*(r4_c(Ni  )/1e12)/20.d0

    gpot_cmb = -(rho0_oc*(r2_c(N)/1e6)/6.d0      + rho1_oc*(r3_c(N)/1e9)/12.d0 + &
                 rho2_oc*(r4_c(N)/1e12)/20.d0    + rho3_oc*(r5_c(N)/1e15)/30.d0)! - &
!               gpot_icb + gpot_icbi)

    rho_psicmb_cmb = gpot_cmb*(rho0_oc*(r3_c(N)/1e9)/3.d0   + rho1_oc*(r4_c(N)/1e12)/4.d0  + &
                               rho2_oc*(r5_c(N)/1e15)/5.d0  + rho3_oc*(r6_c(N)/1e18)/6.d0)
    rho_psicmb_icb = gpot_cmb*(rho0_oc*(r3_c(Ni)/1e9)/3.d0  + rho1_oc*(r4_c(Ni)/1e12)/4.d0 + &
                               rho2_oc*(r5_c(Ni)/1e15)/5.d0 + rho3_oc*(r6_c(Ni)/1e18)/6.d0)

    five_coeff   = rho0_oc*rho0_oc/30.d0
    six_coeff    = rho0_oc*rho1_oc/24.d0
    seven_coeff  = rho1_oc*rho1_oc/84.d0  + rho0_oc*rho2_oc*13.d0/420.d0
    eight_coeff  = rho1_oc*rho2_oc/60.d0  + rho0_oc*rho3_oc/40.d0
    nine_coeff   = rho2_oc*rho2_oc/180.d0 + 7.d0*rho1_oc*rho3_oc/540.d0
    ten_coeff    = rho2_oc*rho3_oc/120.d0
    eleven_coeff = rho3_oc*rho3_oc/330.d0

    rho_psi_cmb = five_coeff*(r5_c(N)/1e15) + six_coeff*(r6_c(N)/1e18)  + seven_coeff*(r7_c(N)/1e21) + &
                  eight_coeff*(r8_c(N)/1e24)+ nine_coeff*(r9_c(N)/1e27) + ten_coeff*(r10_c(N)/1e30)  + &
                  eleven_coeff*(r11_c(N)/1e33)

!    five_coeff   = rho0_ic*rho0_oc/30.d0
!    six_coeff    = rho0_ic*rho1_oc/36.d0
!    seven_coeff  = rho0_ic*rho2_oc/42.d0  + rho2_ic*rho0_oc/140.d0
!    eight_coeff  = rho0_ic*rho3_oc/48.d0  + rho2_ic*rho1_oc/160.d0
!    nine_coeff   = rho2_ic*rho2_oc/180.d0
!    ten_coeff    = rho2_ic*rho3_oc/200.d0

    rho_psi_icb = five_coeff*(r5_c(Ni)/1e15) + six_coeff*(r6_c(Ni)/1e18)  + seven_coeff*(r7_c(Ni)/1e21) + &
                  eight_coeff*(r8_c(Ni)/1e24)+ nine_coeff*(r9_c(Ni)/1e27) + ten_coeff*(r10_c(Ni)/1e30)  + &
                  eleven_coeff*(r11_c(Ni)/1e33)

    int_rhopsi_dV = 16.d0 * pi * pi * G * (rho_psi_cmb - rho_psi_icb + rho_psicmb_cmb - rho_psicmb_icb) *1e15

    Qg = ( Cp + Cc*Cr )  * alphac * (int_rhopsi_dV - M_oc*psir_c(Ni))

    Eg = Qg/Tc

  end subroutine gravitational_poly_core

  subroutine adiabatic_heating(N, Ni, L, Tc, Cr, Cc, alphac, QP, QPL, EP, EPL)

     integer, intent(in)           :: N, Ni
     double precision, intent(in)  :: Tc, Cr, Cc, alphac, L(0:N)
     double precision, intent(out) :: QP, QPL, EP, EPL
     double precision              :: fac, roc(N-Ni)
     double precision              :: integrand_rhor_c2(N-Ni), integrand_rho_rm2(N-Ni), integrand(N-Ni), P_T(N-Ni)
     integer                       :: i

     fac        = -8.d0 * pi * G * alphac * Cr * Cc

     integrand_rhor_c2 = 0.d0
     do i = Ni+1, N
       integrand_rhor_c2(i-Ni) = rhor_c(i)*r_c(i)*r_c(i)
       roc(i-Ni)             = r_c(i)
     enddo
     call splat(integrand_rhor_c2,roc,(N-Ni))

     integrand_rho_rm2 = 0.d0
     do i = Ni+1, N
       integrand_rho_rm2(i-Ni) = integrand_rhor_c2(i-Ni) * rhor_c(i) /r2_c(i)
     enddo
     call splat(integrand_rho_rm2,roc,(N-Ni))

     integrand = 0.d0
     do i = Ni+1, N
       integrand(i-Ni)         = integrand_rho_rm2(i-Ni) - integrand_rho_rm2(N-Ni)
     enddo

     P_T = 0.d0
     P_T = fac * integrand

     integrand = 0.d0
     do i = Ni+1, N
       integrand(i-Ni) = alphaT_c * P_T(i-Ni) * Tar_c(i)* r2_c(i)
     enddo
     call splat(integrand,roc,(N-Ni))

     QP = 4.d0 * pi * integrand(N-Ni)

     write(*,'(A16,E16.6)') 'Qpt=',QP

     integrand = 0.d0
     do i = Ni+1, N
       integrand(i-Ni) = alphaT_c * P_T(i-Ni) * r2_c(i)
     enddo
     call splat(integrand,roc,(N-Ni))

     EP  = QP/Tc - (4.d0 * pi * integrand(N-Ni))

     write(*,'(A16,E16.6)') 'Ept=',EP

     QPL = -4.d0 * pi * r2_c(Ni+1) * rhor_c(Ni+1) * dTmdP_c(Ni+1) * P_T(1) * Cr * L(Ni+1)

     write(*,'(A16,E16.6)') 'QPLt=',QPL

     EPL = QPL * (1.d0/Tc - 1.d0/Tar_c(Ni+1))
     write(*,'(A16,E16.6)') 'Ept=',EPL

  end subroutine adiabatic_heating

  subroutine heatofreaction_num_correct_core(N,Ni,dridt,Cc,Eh)

    integer         , intent(in)  :: N, Ni
    double precision, intent(in)  :: dridt,Cc
    double precision, intent(out) :: Eh
    double precision              :: integrand(N+1), dmudT(0:N)
    integer                       :: i

    dmudT = dmudT0 + dmudT1*r_c/1e3

    dmudT = dmudT * Ev * Na * 1000.d0 / (AO*1.d0)

    integrand = 0.d0
    do i = 0, N
      integrand(i+1) = r_c(i)*r_c(i)*rhor_c(i)*dmudT(i)
    enddo
    call splat(integrand,r_c,N+1)

    Eh = -4.d0 * pi * Cc * dridt * ( integrand(N) - integrand(Ni+1) )

  end subroutine heatofreaction_num_correct_core

! ASSUMES - alpha_c and alpha_D are constants.
!           core is well-mixed and thermodiffusion is negligible (grad(c) & K_T = 0 in eq 15 of Gub etal 14)
!           IS THERE AN ERROR IN EQN 32 OF GUBBINS ET AL 2004?
  subroutine baro_entropy_core(Nic, cbarO_oc, cbarS_oc, cbarSi_oc, Ealpha)

    integer         , intent(in):: Nic
    double precision, intent(in):: cbarO_oc, cbarS_oc, cbarSi_oc
    double precision            :: alphaDO, alphaDS, alphaDSi
    double precision            :: dmubar_dcbarO_oc, dmubar_dcbarS_oc, dmubar_dcbarSi_oc
    double precision            :: dmudcO       , dmudcS       , dmudcSi
    double precision            :: Tbar, rhobar, Abar
    double precision            :: integrand(N+1)
    double precision            :: EalphaO, EalphaS, EalphaSi, Ealpha
    integer                     :: i

    Abar = cbarO_oc*AO + cbarS_oc*AS + cbarSi_oc*ASi + (1.d0-cbarO_oc-cbarS_oc-cbarSi_oc)*AFe

    Tbar   = (Tar_c(N) + Tar_c(0)) / 2
    rhobar = (rhor_c(N) + rhor_c(0)) / 2

    dmubar_dcbarO_oc  = kb*Tbar/cbarO_oc  + lambdaO_oc
    dmubar_dcbarS_oc  = kb*Tbar/cbarS_oc  + lambdaS_oc
    dmubar_dcbarSi_oc = kb*Tbar/cbarSi_oc + lambdaSi_oc

    dmudcO  = dmubar_dcbarO_oc  * (Abar/AO)  * Ev * Na * (1000.d0/AO)
    dmudcS  = dmubar_dcbarS_oc  * (Abar/AS)  * Ev * Na * (1000.d0/AS)
    dmudcSi = dmubar_dcbarSi_oc * (Abar/ASi) * Ev * Na * (1000.d0/ASi)

    alphaDO  = rhobar * DO  / dmudcO
    alphaDS  = rhobar * DS  / dmudcS
    alphaDSi = rhobar * DSi / dmudcSi

    integrand = 0.d0
    do i = 0, N
      integrand(i+1) = r_c(i)*r_c(i)*gr_c(i)*gr_c(i)*(1.d0/Tar_c(i))
    enddo
    call splat(integrand,r_c,N+1)

    EalphaO  = 4.d0* pi * alphac_cO *alphac_cO *alphaDO  * (integrand(N)-integrand(Nic+1))
    EalphaS  = 4.d0* pi * alphac_cS *alphac_cS *alphaDS  * (integrand(N)-integrand(Nic+1))
    EalphaSi = 4.d0* pi * alphac_cSi*alphac_cSi*alphaDSi * (integrand(N)-integrand(Nic+1))

    Ealpha   = EalphaO + EalphaS + EalphaSi

    write(*,'(A16,2E16.6)') 'EalphaO=',EalphaO, alphaDO
    write(*,'(A16,2E16.6)') 'EalphaS=',EalphaS, alphaDS
    write(*,'(A16,2E16.6)') 'EalphaSi=',EalphaSi, alphaDSi

  end subroutine baro_entropy_core

!******************************************************************************
! conductivity_poly: Thermal conductivity using equation (17)
!       of Davies (2014, submitted to PEPI)
!******************************************************************************
  subroutine conductivity_poly_core()

    integer :: i

    do i = 0, N
      kr_c(i) =  k2_oc*r_c(i)*r_c(i)/1e6 + k1_oc*r_c(i)/1e3 + k0_oc
    enddo

    write(*,'(A16,E16.6)') 'k(ro)=',kr_c(N)

  end subroutine conductivity_poly_core

  subroutine conductivity_T_core()

    integer :: i
    double precision :: rho(0:N), Tc(0:N)
    double precision, parameter :: Lo = 2.44e-8

    !Tc   = Tar_c - 273.15 !Kelvin->Celsius
    !rho  = (k0_oc*Tc  + k1_oc) * 1e-8

    !kr_c =  Lo * Tar_c / rho
    kr_c = k0_oc  + k1_oc*Tar_c
    kr_c = kr_c - (pr_c-10e9)*kr_c*0.1/30e9

    write(*,'(A16,2E16.6)') 'k(ro)=', kr_c(N), kr_c(0)

  end subroutine conductivity_T_core

  subroutine cond_entropy_poly_core(N,Ea)

    integer, intent(in)           :: N
    double precision, intent(out) :: Ea
    double precision              :: integrand(N)
    integer                       :: i

    integrand = 0.d0

    do i = 1,N
        integrand(i)=kr_c(i)*(r_c(i)*dTadr_c(i)/Tar_c(i))**2
    enddo
    call splat(integrand,r_c,N)

    Ea=4.d0*pi*integrand(N)

  end subroutine cond_entropy_poly_core

! Find the position where dTdr-dTadr_c changes sign
  subroutine layer_depth2(dTplusdr, dcplusdr, rs)

    double precision, intent(in)  :: dTplusdr, dcplusdr
    double precision, intent(out) :: rs
    double precision              :: diff, rmrs
    integer :: i, ns

    write(*,*) '***************************'
    write(*,*) '*****LAYER DEPTH CALC2*****'
    write(*,'(A16,E16.6)') 'dTdr = ', dTplusdr
    write(*,'(A16,E16.6)') 'dcdr = ', dcplusdr

    rmrs = 1d10
    diff = 0.d0

    do i = 0, N
       diff = dTplusdr - dTadr_c(i)
       if(dabs(diff) <  dabs(rmrs) ) then
         ns   = i
         rmrs = diff
       endif
    enddo

    rs = r_c(ns)

    write(*,'(A16,F16.6,I16)') 'rs   = ', rs, ns
    write(*,*) '***************************'

  end subroutine layer_depth2


! Finds the temperature at (near) the base of the stable layer
  subroutine find(N,rs,A,ns,Tc)

    implicit none

    integer         , intent(in)  :: N
    double precision, intent(in)  :: rs, A(N)
    double precision, intent(out) :: Tc
    double precision              :: rmrs, rdiff
    integer                       :: i, ns

    rmrs  = 1d10
    rdiff = 0.d0

    do i = 0, N
       rdiff = r_c(i) - rs
       if(dabs(rdiff) <  dabs(rmrs) ) then
         ns   = i
         rmrs = rdiff
       endif
    enddo

    if(ns > N) ns = N
    ns = ns+1

    Tc = A(ns)

  end subroutine find

  subroutine open_files_core(char)

    character (len=4)   :: char
    character (len=100) :: energy_file, entropy_file, diag_file, pro_core_file, icb_file, conc_file, snow_file, pro_snow_file

    energy_file  = out_file(1:out_file_len)//"_"//char//'_energy'
    entropy_file = out_file(1:out_file_len)//"_"//char//'_entropy'
    diag_file    = out_file(1:out_file_len)//"_"//char//'_diagnostics'
    pro_core_file= out_file(1:out_file_len)//"_"//char//'_profiles_core'
    icb_file     = out_file(1:out_file_len)//"_"//char//'_icb'
    conc_file    = out_file(1:out_file_len)//"_"//char//'_conc'
    snow_file    = out_file(1:out_file_len)//"_"//char//'_snow'
    pro_snow_file= out_file(1:out_file_len)//"_"//char//'_profiles_snow'

    open(10,file=energy_file)
    open(11,file=entropy_file)
    open(12,file=diag_file)
    open(13,file=pro_core_file)
    open(14,file=icb_file)
    open(15,file=conc_file)
    open(16,file=snow_file)
    open(17,file=pro_snow_file)
    write(10,'(12A16)') 'Time (Myr)','Qs','Qg','Ql','Qr','Qcmb','Qa'              ,'QL_rs','QL_snow','QL_melt','Qg_snow','Qg_melt'
    write(11,'(14A16)') 'Time (Myr)','Es','Eg','El','Er','EJ'  ,'Ea','Eh','Ealpha','EL_rs','EL_snow','EL_melt','Eg_snow','Eg_melt'
    write(12,'(9A16)')  'Time (Myr)', 'Tc', 'dTcdt', 'ri', 'rs', 'dridt', 'M_c', 'M_oc','M_snow'
    write(13,'(9A16)')  'r', 'T', 'dTadr_c', 'Tm', 'p', 'rho', 'g', 'psi', 'k'
    write(14,'(11A16)') 'r', 'dTadP(ri)','dTmdP(ri)','rho(ri)','T(ri)','g(ri)','p(ri)','Cr','dTadr(ri)','dTmdr(ri)','T_cen'
    write(15,'(14A16)') 'Time (Myr)','ri','clO','clS','clSi','csO','csS','csSi','QgO','QgS','QgSi','dTmO', 'dTmS', 'dTmSi'
    write(16,'(10A16)') 'Time (Myr)','rd', 'drsdt', 'Cp', 'Cr', 'CcS', 'Ql', 'Qg', 'El', 'Eg'
    write(17,'(6A16)')  'r', 'xil', 'phi', 'dcldT_s', 'Cl_snow', 'L'

  end subroutine open_files_core

  subroutine close_files_core()

    close(10)
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(17)

  end subroutine close_files_core

  subroutine write_profiles_core()

    implicit none

    integer :: j

    do j = 0, N
      write(13,'(9E16.8)') r_c(j), Tar_c(j), dTadr_c(j), Tm_c(j), pr_c(j), rhor_c(j), gr_c(j), psir_c(j), kr_c(j)
    enddo
    write(13,*)
    write(13,*)

  end subroutine write_profiles_core

  subroutine write_profiles_snow(N, Nsl, time, clS, cliq_s, dcldT_s, phi_s, LH_c, Cl_snow)

    integer, intent(in)          :: N, Nsl
    double precision, intent(in) :: time, clS, cliq_s(0:N), dcldT_s(0:N), phi_s(0:N), LH_c(0:N), Cl_snow(0:N)
    integer                      :: i 
 
    write(17,*) time
    do i = 0, Nsl
      write(17,'(6E16.8)') r_c(i), clS+cliq_s(i), phi_s(i), dcldT_s(i), Cl_snow(i), LH_c(i)
    enddo
    do i = Nsl+1, N
      write(17,'(6E16.8)') r_c(i), cliq_s(i)    , phi_s(i), dcldT_s(i), Cl_snow(i), LH_c(i)
    enddo
    write(17,*)
    write(17,*)

  end  subroutine write_profiles_snow

end module core_model
