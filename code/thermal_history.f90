  program thermal_history

    use core_evolution
    use mantle_evolution

    implicit none

    double precision :: Q, Tc, dTcdt, ri, time, EJ
    double precision :: cbarO_oc, cbarS_oc, cbarSi_oc, cO, cS, cSi
    integer          :: forback, ncla
    character(len=4) :: forbackc
    character(len=200) :: arg

    ncla = command_argument_count()

    call getarg(1,arg)
    arg = trim(arg)
    call read_input_core(arg,Q,Tc,ri,EJ,cbarO_oc, cbarS_oc, cbarSi_oc) ! Read the input file
    call rpts_core(ri,rs)                                          ! Create grid so Tm and Ta can be obtained

    if(ncla==2) then 
      call getarg(2,arg)
      arg = trim(arg)
      call read_input_mantle(arg)   
    elseif(ncla>2) then
      stop ">2 CLAs"
    endif

    !Tm   = Tc
    Ms1  = 0.d0; Ms2  = 0.d0; Mc1  = 0.d0; Mc2  = 0.d0
    cMs1 = 0.d0; cMs2 = 0.d0; cMc1 = 0.d0; cMc2 = 0.d0
    ri0 = ri  ; Nic = Ni

    iteration = 1
    forback   = 1                                                  ! Tells code to go backwards in time (forwards =1)
    forbackc  = "forw"
    call open_files_core(forbackc)                                 ! Open output files
    call open_files_mantle(forbackc)                               ! Open output files

    do i = 0, nt

      call get_time(i,time,dt,forback)
      write(*,*) '----------------------', i

      Q = cmbflux(i)
      if(sol==0) then
        time = time_file(nt-i)
        dt   = dt_file(nt-i)
        Q    = cmbflux(nt-i)
      elseif(sol==4) then 
        Q    = cmbflux(0) + cmbflux(1)*dexp(-time*1e-3/cmbflux(2))
      elseif(sol==6) then 
        call mantle_evol_calc(Tc,rc,time,forback,Tm,Q)
      elseif(sol==7) then
        call mantle_evol_calc_breuer(Tc,rc,time,forback,Tm,Q)
      endif
      if(i==0) iteration = 1

      call mole2massconc_core(cbarO_oc, cbarS_oc, cbarSi_oc, cO, cS, cSi)
      call core_evol_calc(Q,Tc,time,dTcdt,ri,EJ,cO,cS,cSi,forback)   !   arg1=forward(1) or backward(2)
      call mass2moleconc_core(cO, cS, cSi, cbarO_oc, cbarS_oc, cbarSi_oc)

      ri0 = r_c(Nic)
  
      if(forback < 0 ) goto 100
      if(forback == 3) goto 100
    enddo

100 continue

    if(forback==3) write(*,*) 'Snow zone occupies whole core - Terminating'

    call close_files_core()

    write(*,*) '---END OF FORWARD  CALC---'
    write(*,'(A16,F16.6)') 'Tcen = ', T_cen
    write(*,'(A16,F16.6)') 'Tc   = ', Tc
    write(*,'(A16,F16.6)') 'ri   = ', ri
    write(*,'(A16,E16.6)') 'Qcmb = ', Q

    deallocate(cmbflux)

  end program thermal_history
