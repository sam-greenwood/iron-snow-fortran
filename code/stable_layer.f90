  module stable_layer

    use parameters_core

    implicit none
    save

    integer, parameter  :: num_file=100, an1_file=111, an2_file=222
    integer, parameter  :: Ns_tot=N, Ns_it=20 !# pts added per iteration
    
    double precision    :: Ta_s(Ns_tot), dTadr_s(Ns_tot), x_s(Ns_tot)
    double precision    :: courant

  contains

!-------------------------------------------------------------
! Driver program
!
! Inputs: xmin: min radius of stable layer (metres)
!         xmax: max radius of stable layer (metres)
!         Diff: Diffusivity (m^2/s)
!         bval: value on bottom boundary
!         tval: value on top boundary
!         time0: initial time
!         ttot: Total time to integrate
!         scalar: temperature or concenration
!         SOL: on input-IC; on output-soln
!
!   NOTE: prescribing Courant means that different resolutions will give
!     different answers even when 2 runs have taken the same # steps because
!     the step length (dt) is different
!     This also happens when you change D; the result is changing dt for
!     fixed dx (as specified by the grid) and courant
!-------------------------------------------------------------

 subroutine stable_layer_calc(Nstrat,xmin,xmax,Diff,bval,bbc,tval,tbc,S,time0,ttot,scalar,SOL)

    integer         , intent(in)    :: Nstrat
    double precision, intent(in)    :: xmin,xmax,Diff,bval,tval,time0,ttot,S(Nstrat)
    double precision, intent(inout) :: SOL(Nstrat)
    double precision                :: DL(Nstrat-1),D(Nstrat),DU(Nstrat-1),RHS(Nstrat),x(Nstrat), dx
    double precision                :: dt, time,Src(Nstrat), SOLGUB(Nstrat)
    integer                         :: i, INFO, Nt, j,bbc,tbc
    character (len=4)               :: scalar

    write(*,*) '*****STABLE LAYER CALC*****'
    write(*,'(A16,I16)') 'N     = ', Nstrat
    write(*,'(A16,E16.6)') 'xmin  = ', xmin 
    write(*,'(A16,E16.6)') 'xmax  = ', xmax
    write(*,'(A16,E16.6)') 'D     = ', Diff
    Write(*,'(A16,E16.6)') 'Smax  = ', maxval(S)
    write(*,'(A16,E16.6,I16)') 'xmax BC = ', tval,tbc
    write(*,'(A16,E16.6,I16)') 'xmin BC = ', bval,bbc
    write(*,'(A16,E16.6)') 'time  = ', time0

    if( dabs(xmin-xmax) < 1d-3) then 
      write(*,*) 'xmin = xmax! Returning...'
      return
    endif

    ! Set up a grid for the core
    call grid(Nstrat, xmin, xmax, dx, x)

!   Results file
    call open_strat_files(scalar,time0,num_file)
    write(num_file,'(2A16)') 'x', scalar

    courant = 0.4d0
    dt      = courant*dx*dx / Diff                  !Courant tells you timestep
!    dt = ttot/1000.d0
    Nt      = ttot/dt                               !# timepts needed for 1DT

    if(Nt==0) then 
      write(*,*) 'Taking 1 step even though Nt=0...'
      Nt=1
    endif

    write(*, '(A16,2E16.6)') 'Delta t,x        = ', dt, dx
    write(*, '(A16,I16)')   '# tpts for 1DT = ', Nt
    write(*,*) '***************************'

    call write( Nstrat,i,time,num_file,x,SOL)                          !Write solution

    do i = 1, Nt

      Src  = 2d0*dt*S
      RHS  = 0.d0
      time = time0 + real(i)*dt

      call lhs_construct(Nstrat,courant,DL,D,DU)         !Construct matrix elements
      call lhs_bc(Nstrat,bbc,1,DL,D,DU)                  !LHS Bottom BC (arg3=1) (arg2=1 for Arfken)
      call lhs_bc(Nstrat,tbc,2,DL,D,DU)                  !LHS Top    BC (arg3=2) (arg2=1 for Arfken)
      call rhs_construct(Nstrat,courant,SOL,Src,RHS)     !Fill the rest of the RHS
      call setbc (Nstrat,bbc,1,bval,courant,dx,SOL,RHS)
      call setbc (Nstrat,tbc,2,tval,courant,dx,SOL,RHS)

      call dgtsv(Nstrat, 1, DL, D, DU, RHS, Nstrat, INFO )  !Get T(t_i+1)
      if(INFO .ne. 0) then
        write(*,*) "Call:       DGBTRS"
        write(*,*) "INFO      = ", INFO
        stop
      endif

!      call analytical_soln1(Nstrat,x,time,Diff,tval,SOLGUB)
      call write( Nstrat,i,time,num_file,x,RHS)                    !Write solution
!      call write( Nstrat,i,time,567  ,x,SOLGUB)                    !Write solution

      SOL = RHS

    enddo

  end subroutine stable_layer_calc

!-------------------------------------------------------------
! Set up a grid of evenly-spaced points from xmin->xmax
! xmin/xmax in km
!-------------------------------------------------------------

  subroutine grid(Nstrat, xmin, xmax, dx, x)

    implicit none

    integer, intent(in)  :: Nstrat
    integer              :: i
    double precision, intent(in)  :: xmin, xmax
    double precision, intent(out) :: x(Nstrat), dx


    x  = 0.d0
    dx = (xmax - xmin) / (Nstrat-1)

    x(1) = xmin
    do i = 2, Nstrat
      x(i) = x(i-1) + dx
    enddo

  end subroutine grid

  subroutine lhs_construct(Nstrat,coeff,DL,D,DU)

    implicit none

    integer         , intent(in)  :: Nstrat
    double precision, intent(in)  :: coeff
    double precision, intent(out) :: DL(Nstrat-1), D(Nstrat), DU(Nstrat-1)

    DL = 0.d0; D = 0.d0; DU = 0.d0
    DL = -coeff
    D  = 2.d0 + 2*coeff
    DU = -coeff

  end subroutine lhs_construct

!-------------------------------------------------------------
! Boundary conditions
! bc == 1 => Dirichlet; ==2 => Neumann
! tb == 1 => bottom; ==2 => top
!-------------------------------------------------------------

  subroutine lhs_bc(Nstrat,bc,tb,DL,D,DU)

    implicit none

    integer         , intent(in)  :: Nstrat, bc, tb
    double precision              :: DL(Nstrat-1), D(Nstrat), DU(Nstrat-1)

    if(bc==1 .and. tb==1) then 
      D(1)  = 1.d0
      DU(1) = 0.d0
    elseif(bc==1 .and. tb==2) then
      D(Nstrat)    = 1.d0
      DL(Nstrat-1) = 0.d0
    elseif(bc==2 .and. tb==1) then !D is correct from lhs_construct
      DU(1)    = 2.d0*DU(1) 
    elseif(bc==2 .and. tb==2) then
      DL(Nstrat-1) = 2.d0*DL(Nstrat-1)
    else 
      stop "LHS BC: BC's must be Dirichlet or Neumann!"
    endif

  end subroutine lhs_bc

!-------------------------------------------------------------
! Set the top and bottom rows of the RHS
! NOTE - rhs_construct and setbc both use T^n. Need to ensure
! they are not using modified values, i.e. rhs_construct mods 
! The whole RHS vector, but setbc needs to use the original 
! T^n's
!-------------------------------------------------------------

  subroutine setbc(Nstrat, bc, tb, ans, coeff, dx, SOL, RHS)

    implicit none

    integer         , intent(in)    :: Nstrat, bc, tb
    double precision, intent(in)    :: ans, coeff, dx
    double precision, intent(inout) :: RHS(Nstrat), SOL(Nstrat)

    if(bc==1 .and. tb==1) then     !Bottom, Fixed T
      RHS(1)  = ans
    elseif(bc==1 .and. tb==2) then !Top, Fixed T
      RHS(Nstrat) = ans
    elseif(bc==2 .and. tb==1) then 
      RHS(1)      = (2.d0 - 2.d0*coeff)*SOL(1)      + 2.d0*coeff*SOL(2)        - 4.d0*coeff*ans*dx
    elseif(bc==2 .and. tb==2) then
      RHS(Nstrat) = (2.d0 - 2.d0*coeff)*SOL(Nstrat) + 2.d0*coeff*SOL(Nstrat-1) + 4.d0*coeff*ans*dx
    else
      stop "RHS BC: BC's must be Dirichlet or Neumann!"
    endif

  end subroutine setbc

!-------------------------------------------------------------
! RHS for Crank-Nicholson scheme
!-------------------------------------------------------------

  subroutine rhs_construct(Nstrat,coeff,SOL,S,RHS)

    implicit none

    integer                         :: i
    integer, intent(in)             :: Nstrat
    double precision, intent(in)    :: coeff, SOL(Nstrat), S(Nstrat)
    double precision, intent(out  ) :: RHS(Nstrat)

    do i = 2, Nstrat-1
      RHS(i) = coeff*SOL(i-1) + (2.d0 - 2.d0*coeff)*SOL(i) + coeff*SOL(i+1) + S(i)
    enddo

  end subroutine rhs_construct

!-------------------------------------------------------------
! Open files
!-------------------------------------------------------------

  subroutine open_strat_files(char,time,unit)

    implicit none

    integer, intent(in)             :: unit
    double precision, intent(in)    :: time
    character (len=4)               :: char
    character (len=5)               :: cNstrat
    character (len=9)               :: ct
    character (len=100)             :: dir

    write (ct, '(E9.3)'), time

    dir   = char//"_t="//ct
    open(unit, file=dir)

  end subroutine open_strat_files

!-------------------------------------------------------------
! Close files
!-------------------------------------------------------------

  subroutine close_files(unit)

    implicit none

    integer, intent(in)             :: unit

    close(unit)

  end subroutine close_files

!-------------------------------------------------------------
! Write data to file
!-------------------------------------------------------------

  subroutine write(Nstrat,it,time,unit,x,T)

    implicit none

    integer          :: Nstrat, unit, it
    double precision :: x(Nstrat), T(Nstrat), time
    integer          :: i

    write(unit,'(A12,E16.6,I12)') 'time = ', time, it
    do i = 1, Nstrat
      if(T(i) < 1d-20) T(i) = 0.d0
      write(unit,*) x(i), T(i)
    enddo
    write(unit,*)
    write(unit,*)

  end subroutine write

  subroutine topbc(D, ans)

    implicit none

    double precision            :: D, ans
    double precision, parameter :: alphac=1.1d0, alphad=0.62d-12, g=10.68d0, rho=9903.49d0

    ans = (alphac*alphad*g)/(rho*D)

  end subroutine topbc


! Use for anaytical soln 3
  subroutine write2(Nstrat,it,time,unit,x,T)

    implicit none

    integer          :: Nstrat, unit, it
    double precision :: x(Nstrat), T(Nstrat), time
    integer          :: i

    write(unit,'(A12,E16.6,I12)') 'time = ', time, it
    do i = 1, Nstrat
      if(T(i) < 1d-20) T(i) = 0.d0
      write(unit,*) x(Nstrat-i+1), T(i)
    enddo
    write(unit,*)
    write(unit,*)

  end subroutine write2

!-------------------------------------------------------------
! This subroutine implements the analytical soln of Gubbins &
! Davies (PEPI, 2013)
! At t = 0, T=0 for ri=3400 < r < ro=3480 (km)
!           T=0 for r=ri
!           dT/dr=-b for r=ro
! Here we will use b=1 for simplicity
!-------------------------------------------------------------

  subroutine analytical_soln1(Nstrat,x,t,D,b,Ta)

    implicit none

    integer          :: Nstrat, i
    double precision :: x(Nstrat), Ta(Nstrat), tmp(Nstrat), t, D, b
    double precision :: h, t1, t2, zeta, dx, dy, y
    double precision, parameter :: pi  = 3.1415926535897931d0

    h    = dsqrt(D*t); Ta = 0.d0; y = 0.d0

    do i = 1, Nstrat
       dy  = x(i) - x(i-1)
       if(i==1) dy = 0.d0
       y   = y + dy

       zeta = y/(2.d0*h)

       t1   = dexp(-(zeta**2)) / dsqrt(pi)
       t2   = zeta*derfc(zeta)

       Ta(i) = 2.d0*b*h*(t1-t2)
    enddo

  end subroutine analytical_soln1


!-------------------------------------------------------------
! This subroutine implements the analytical soln on pg. 612 of
! Arfken (eqn 9.223).
! At t = 0, T=1 for -1 < x < 1
!           T=0 for x=-1, 1
!-------------------------------------------------------------

  subroutine analytical_soln2(Nstrat,x,t,D,Ta)

    implicit none

    integer          :: Nstrat, i, m
    double precision :: x(Nstrat), Ta(Nstrat), t, D, Ds
    double precision :: expfac, tmp1, f1, f2, f3
    double precision, parameter :: pi  = 3.1415926535897931d0

    Ds = dsqrt(D)

    Ta = 0.d0
    do i = 1, Nstrat
      f1 = 0.d0; f2 = 0.d0; f3 = 0.d0
      do m = 0, 1000
         tmp1 = 2.d0*m + 1.d0
         expfac = ((tmp1) * pi * Ds) / 2.d0
         f1 = ((-1.d0)**m) / tmp1
         f2 = dcos(tmp1 * (pi*x(i)/2.d0))
         f3 = dexp(-t*expfac**2)

         Ta(i) = Ta(i) + (f1*f2*f3)
      enddo
    enddo

    Ta = Ta * (4.d0/pi)

  end subroutine analytical_soln2

  subroutine analytical_soln3(Nstrat,x,t,D,b,Ta)

    implicit none

    integer          :: Nstrat, i
    double precision :: x(Nstrat), Ta(Nstrat), tmp(Nstrat), t, D, b
    double precision :: h, t1, t2, zeta, dx, dy, y
    double precision, parameter :: pi  = 3.1415926535897931d0

    h    = dsqrt(D*t); Ta = 0.d0; y = 0.d0

    do i = 1, Nstrat
       dy  = x(i) - x(i-1)
       if(i==1) dy = 0.d0
       y   = y + dy

       zeta = y/(2.d0*h)

       Ta(i) = 2.d0*b*h*erf(zeta)
    enddo

  end subroutine analytical_soln3

! Solution for Neumann conditions at both ends and an initial T of
! 9-3*cos(pi*x/4) - 6*cos(2*pi*x) on a bar of length 8. 

  subroutine analytical_soln4(Nstrat,x,t,D,Ta)

    integer          :: Nstrat, i
    double precision :: x(Nstrat), Ta(Nstrat), term1(Nstrat), term2(Nstrat), term3(Nstrat), t, D
    double precision, parameter :: pi  = 3.1415926535897931d0

    term1 = 9.d0
    term2 = D   * dexp(-(D * 4.d0   * pi**2 * t) / 64.d0) * dcos(pi*x/4.d0)
    term3 = 2*D * dexp(-(D * 256.d0 * pi**2 * t) / 64.d0) * dcos(2.d0*pi*x)

    Ta = term1 - term2 - term3

  end subroutine analytical_soln4

! CJ pg 79 eqn 2 with S = A0/K
  subroutine analytical_soln5(Nstrat,x,t,D,T0,S,Ta)

    integer          :: Nstrat, i
    double precision :: x(Nstrat), Ta(Nstrat), T0, S, t, D, h, y, dy
    double precision :: term1, term2, term1_fac, term2_fac, zeta
    double precision, parameter :: pi  = 3.1415926535897931d0

    Ta = 0.d0; y = 0.d0

    h         = dsqrt(D*t); Ta = 0.d0; y = 0.d0
    term1 = 0.d0; term2 = 0.d0

    do i = 1, Nstrat
       dy  = x(i) - x(i-1)
       if(i==1) dy = 0.d0
       y   = y + dy

       zeta = y/(2.d0*h)

       term1_fac = (T0 + D*t*S + S*y**2/2.d0)
       term2_fac = S*y*dsqrt(D*t/pi)

       term1 = term1_fac * derf(zeta)
       term2 = term2_fac * dexp(-(zeta**2))

       Ta(i) = term1 + term2 - S*y**2/2.d0

    enddo

  end subroutine analytical_soln5

  subroutine adiabat_strat_nimmo(N, x, Tc, Ta, dTadr, Qa)

      double precision, parameter  :: Da=0.5d7

      integer         , intent(in) :: N
      double precision, intent(in) :: Tc, x(N)
      double precision, intent(out):: Ta(N), dTadr(N)
      double precision :: Qa

      Ta    = Tc*dexp(-x**2/Da**2)
      dTadr = -2.d0*x *Tc*dexp(-x **2/Da**2)/Da**2

    end subroutine adiabat_strat_nimmo


end module stable_layer
