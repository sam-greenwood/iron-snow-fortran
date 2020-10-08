  module core_evolution

    use core_model
    use parameters

    implicit none

    ! Note that pressure doesn't change in time, but this doesn't matter
    ! as it only affects the melting curve, which is assumed const
    ! t indicates 'tilde' i.e. not yet multiplied by dTcdt

    ! The pressure could, in principle, be evolved in time because it is
    ! Inside the loops rather than precomputed. This is because we need all
    ! the profiles on the same grid for computing the ICB radius in ricalc

    integer          :: Nic, Nsl, i, j, z
    double precision :: Qa, Qst, Qr, Qlt, Qgt, Qtt, QPt, QPLt
    double precision :: Ealpha, Ea, Est, Er, Elt, Egt, Ett, Eht, EPt, EPLt
    double precision :: Qst_snow, Qlt_snow, Qgt_snow, Qmt_snow, Qgt_liq
    double precision :: Est_snow, Elt_snow, Egt_snow, Emt_snow, Egt_liq
    double precision :: dTplusdr, dcplusdr, dTcmbdr, dTrsdr
    double precision :: drsdt, dridt, Cr_snow, Cr
    double precision :: Ttop, Tbottom, difftb
    double precision :: Ms1, Ms2, Mc1, Mc2, cMs1, cMs2, cMc1, cMc2
    double precision, parameter :: ONES(0:N) = 1d0

  contains

  subroutine core_evol_calc(Q,Tc,time,dTcdt,ri,EJ,clO,clS,clSi,forback)

    double precision :: Q, Tc, time, dTcdt, ri, EJ
    double precision :: QgO, EgO, QgS, EgS, QgSi, EgSi
    double precision :: dMsdt, dMcdt, dcMsdt, dcMcdt
    double precision :: c0_s(0:N), cliq_s(0:N), dcldT_s(0:N), phi_s(0:N), LH_c(0:N), LH_s(0:N)
    double precision :: rm1, rm2, rb
    double precision :: CcO, CcS, CcSi, dTmO(0:N), dTmS(0:N), dTmSi(0:N)
    double precision :: clO, clS, clSi, csO, csS, csSi, CcS_snow, Cp_snow, Cl_snow(0:N), Cl_melt(0:N), xi_melt(0:N), phi_tot
    double precision :: clbarO, clbarS, clbarSi, csbarO, csbarS, csbarSi
    integer          :: forback, j, i

    call mass2moleconc_core(clO, clS, clSi, clbarO, clbarS, clbarSi) !Get mole fractions too
    call adiabat_poly_core(N,Tc,Qa)                                  !Find Ta on this grid
    call pressure_poly_core(Nic,N,Pc)                                !Find P on this grid
    !call conductivity_T_core()                                       !Get conductivity profile
    call conductivity_poly_core()
    call entropy_melting()                                           !Find dS
    call melting_mars(clS)                                       !Find Tm(P)

    if(meltdep == 0.d0) then
      call solid_conc(clbarO , lambdaO_oc , lambdaO_ic , dmuO , csbarO)
      call solid_conc(clbarS , lambdaS_oc , lambdaS_ic , dmuS , csbarS)
      call solid_conc(clbarSi, lambdaSi_oc, lambdaSi_ic, dmuSi, csbarSi)
      call melting_pt_dep(clbarO, csbarO, dTmO)    !Melting pt depression due to O
      call melting_pt_dep(clbarS, csbarS, dTmS)    !Melting pt depression due to S
      call melting_pt_dep(clbarSi, csbarSi, dTmSi) !Melting pt depression due to Si
      call mole2massconc_core(csbarO, csbarS, csbarSi, csO, csS, csSi)
      Tm_c = Tm_c + dTmO + dTmS + dTmSi
    else
      csO=clO; csS=0.0d0*clS; csSi=clSi
    endif

    if( (iteration==0) ) call adiabat_init(Tc,Qa)

    call ricalc(N, ri)                               !Find inner core radius using Tm & Ta on old grid
    call rscalc(N, rs, rm1, rm2, rb)
    write(356,'(6F20.8)') time, rb, rm2, rm1, rs
    if(rs < 1e3) forback = 3

    Nic=Ni                                           !Ni is total number of pts used for IC
    Nsl=N-1                                          !So its not out of range in calls like crfac
    if(ri == 0.d0) then
      Nic=0; csO=0.d0; csS=0.d0; csSi=0.d0
    endif
    if(float(ceiling(rs)) .ne. rc) Nsl=Ni+Ns

    call terminate(rb, rs, forback)
    call rpts_core(ri,rs)                            !Set new grid based on ri

    call mass_poly_core(Nic,Nsl,N)                 ! Mass just for convecting part
    call mass_poly_correct(ri)
    call mass_poly_core(Nic,Nsl,N)
    call density_poly_core(Nic, clO, alphac_cO)
    call gravity_poly_core(Nic)
    call grav_pot_poly_core(Nic,N)
    call adiabat_poly_core(N,Tc,Qa)
    call pressure_poly_core(Nic,N,Pc)
    !call conductivity_T_core()                     ! Get conductivity profile
    call conductivity_poly_core()
    call entropy_melting()
    call melting_mars(clS)
    if(meltdep == 0.d0) then
      call melting_pt_dep(clbarO, csbarO, dTmO)    !Melting pt depression due to O
      call melting_pt_dep(clbarS, csbarS, dTmS)    !Melting pt depression due to S
      call melting_pt_dep(clbarSi, csbarSi, dTmSi) !Melting pt depression due to Si
      Tm_c = Tm_c + dTmO + dTmS + dTmSi
    endif

    call melting_mars_snowzone(Nsl)                !Set Tm = Ta in snow zone
    call lh_coeff_core(N,LH_c)                     !Set L at all radii
    if(float(ceiling(rs)) .ne. rc) then
      call xliq_snow(Nsl, Tar_c, cliq_s, dcldT_s)  !Get xliq and dxliq/dT from Tm
!      call phi_snow(Nsl, Tc, dTcdt, clS, cliq_s, c0_s, phi_s) !Tina's phi estimate
      call liquidusfac(N,Nsl,LH_c,cliq_s,phi_s, Cl_snow)
      dcldT_s = Cl_snow
      phi_s(0:Nsl) = 0.d0
      phi_s(Nsl+1:N) = Cl_snow(Nsl+1:N)*Tar_c(Nsl+1:N)*dTcdt*secinyr*dt/cliq_s(Nsl+1:N)/Tc !phi on new grid with old dT/dt
      call solid_mass_test(N, Nsl, time, phi_s, phi_tot)
      call crfac(Nsl,Tc,rs,Cr_snow)
      call ccfac(Nsl, cliq_s(Nsl+1), clS , CcS_snow) !Equal to zero
      call cpfac(Nsl, cliq_s, dcldT_s, Tc, Cp_snow)
      call latent_core(  Nsl,LH_c,Tc,Cr_snow,Qlt_snow,Elt_snow)
      Qlt_snow = Qlt_snow*phi_s(Nsl+1); Elt_snow = Elt_snow*phi_s(Nsl+1) !Total LH released by freezing at interface
    else 
      cliq_s = clS; dcldT_s=0.d0; phi_s=0.d0; Cl_snow=0.d0; Cr_snow=0.d0; CcS_snow=0.d0; Cp_snow=0.d0
    endif

    ! Inner core  
    call crfac(Nic,Tc,ri,Cr)                                             !dridt
    call ccfac(Nic, clO , csO,  CcO)
    call ccfac(Nic, clS , csS,  CcS)
    call ccfac(Nic, clSi, csSi, CcSi)                                    !Dc/dt -> dTc/dt
    call latent_core(  Nic,LH_c,Tc,Cr,     Qlt,Elt)                      !Latent heat

    !Liquid core
    Qgt = 0.d0; Egt = 0.d0
    call secular_num(N,Tc,Qst,Est)                                              !Qs with standard cp
!    call secular_poly_core(N,Nic,Tc,Qst,Est) !Secular
!    call gravitational_liq(N,Nic,Tc,Cr     ,CcS     ,0.d0   ,alphac_cS,QgS,EgS)!Qg numerically
    call radiogenic_num_core(N,Tc,h,Qr,Er)                                      !Radiogenic
    call gravitational_poly_core(N,Nic,Tc,Cr,CcO, 0.d0,alphac_cO ,QgO,EgO)      !Gravitational energy
    Qgt=QgO; Egt=EgO
    call gravitational_poly_core(N,Nic,Tc,Cr,CcS, 0.d0,alphac_cS ,QgS,EgS)      !Gravitational energy
    Qgt=Qgt+QgS;  Egt=Egt+EgS
    call gravitational_poly_core(N,Nic,Tc,Cr,CcSi,0.d0,alphac_cSi,QgSi,EgSi)    !Gravitational energy
    Qgt=Qgt+QgSi; Egt=Egt+EgSi
    call heatofreaction_num_correct_core(N,Nic,Cr,CcO,Eht)                      !Heat of reaction
    call cond_entropy_poly_core(N,Ea)                                           !Thermal conduction
    if(float(ceiling(rs)) .ne. rc) then
      ! Latent heat released: \int rho *L * d\phi/dt dV = L dMs/dt
      call secular_snow(N,Nsl,Tc,LH_c,cliq_s,Cl_snow,Qst_snow,Est_snow)
      ! Latent heat absorbed: L(rs) dMs/dt; only L is evaluated at the interface 
      LH_s = LH_c(Nsl+1)                                                        !Set LH to base of snow value
      call secular_snow(N,Nsl,Tc,LH_s,cliq_s,Cl_snow,Qmt_snow,Emt_snow)
      call gravitational_remelting(Nsl, 0, Tc, Cr_snow, CcS_snow, Cp_snow, cliq_s(Nsl+1)-clS, alphac_cS, Qgt_liq, Egt_liq) !Changed to be jump in composition-SamG.
      call gravitational_freezing(N,Nsl,Tc,Cr_snow,cliq_s        ,dcldT_s,alphac_cS,Qgt_snow,Egt_snow) !Snow zone
      Qmt_snow = -1.d0*Qmt_snow; Emt_snow = -1.d0*Emt_snow 
    endif
    call baro_entropy_core(Nic, clbarO, clbarS, clbarSi, Ealpha)

    if(ah==1) then                                    !Small QP and QPL terms
      call adiabatic_heating(N, Ni, LH_c, Tc, Cr, CcO, alphac_cO, QPt, QPLt, EPt, EPLt)
    else
      QPt=0.d0; QPLt=0.d0; EPt=0.d0; EPLt=0.d0
    endif
    if(hor .ne. 1) Eht = 0.d0

    call core_en_ent_calc(Q,time,dTcdt,ri,EJ,forback) !Relate EJ and Q
    call write_profiles_core()
    call write_profiles_snow(N, Nsl, time, clS, cliq_s, dcldT_s, phi_s, LH_c, Cl_snow)

    dridt = Cr     *dTcdt                             !IC growth rate
    drsdt = Cr_snow*dTcdt
    dTcmbdr = -Q/(4.d0*pi*rc*rc*kr_c(N))              !dTdr @ CMB
    write(*,'(A16,E16.6)') 'drsdt = ', drsdt

    call le_mass_test(N, Nsl, time,clS, cliq_s)

    call dcdt(CcO      , Cr     , 0.d0,    dTcdt, clO , forback)
    call dcdt(CcS      , Cr     , 0.d0,    dTcdt, clS , forback)
    call dcdt(CcSi     , Cr     , 0.d0,    dTcdt, clSi, forback)
    call dcdt(CcS_snow , Cr_snow, Cp_snow, dTcdt, clS , forback)
    call radiogenic_decay_core(time, h0_c, h)

    write(10,'(12E16.6)') time, Qst*dTcdt, Qgt*dTcdt, Qlt*dTcdt, Qr, &
                          Q, Qa, Qlt_snow*dTcdt,Qst_snow*dTcdt,Qmt_snow*dTcdt,Qgt_snow*dTcdt,Qgt_liq*dTcdt
    write(11,'(14E16.6)') time, Est*dTcdt, Egt*dTcdt, Elt*dTcdt, Er, EJ, Ea, Eht*dTcdt, &
                          Ealpha,Elt_snow*dTcdt,Est_snow*dTcdt,Emt_snow*dTcdt,Egt_snow*dTcdt,Egt_liq*dTcdt
    write(12,'(9E16.8)' ) time, Tc, dTcdt*secingyr,ri, rs, dridt*secingyr, M_c,M_oc, M_s
    write(15,'(14E16.6)') time, ri, clO, clS, clSi, csO, csS, csSi, QgO*dTcdt,QgS*dTcdt, QgSi*dTcdt, &
                          dTmO(Ni+1), dTmS(Ni+1), dTmSi(Ni+1)
    write(16,'(10E16.8)') time,(rc-rs),drsdt*secingyr,Cp_snow,Cr_snow,CcS_snow,Qlt_snow*dTcdt, &
                          Qgt_liq*dTcdt,Elt_snow*dTcdt,Egt_liq*dTcdt

    if(forback==1) Tc  = Tc + dTcdt*secinyr*dt        !New CMB temperature
    if(forback==2) Tc  = Tc - dTcdt*secinyr*dt        !New CMB temperature

    if(Tc < 0.d0) then
      write(*,*) '********************'
      write(*,*) 'STOP: CMB temperature < 0!'
      write(*,*) 'time, Tc = ', time, Tc
      write(*,*) '********************'
      stop
    endif

  end subroutine core_evol_calc

  subroutine core_en_ent_calc(Q,time,dTcdt,ri,EJ,forback)

    double precision,save :: Qlast, Efirst
    double precision      :: time, ri, Q, dTcdt, EJ
    integer               :: forback
    integer, save         :: qq

!     CHK!
      Qtt = Qst + Qgt + Qlt +       QPt + QPLt + Qgt_liq + Qlt_snow + Qst_snow + Qmt_snow + Qgt_snow
      Ett = Est + Egt + Elt + Eht + EPt + EPLt + Egt_liq + Elt_snow + Est_snow + Emt_snow + Egt_snow

      write(*,'(A16,E16.6)') 'Qst=',Qst
      write(*,'(A16,E16.6)') 'Qgt=',Qgt
      write(*,'(A16,E16.6)') 'Qlt=',Qlt
      write(*,'(A16,E16.6)') 'Qst_snow=',Qst_snow
      write(*,'(A16,E16.6)') 'Qmt_snow=',Qmt_snow
      write(*,'(A16,E16.6)') 'Qlt_snow=',Qlt_snow
      write(*,'(A16,E16.6)') 'Qgt_liq =',Qgt_liq
      write(*,'(A16,E16.6)') 'Qgt_snow=',Qgt_snow
      write(*,'(A16,E16.6)') 'Qtt=',Qtt
      write(*,'(A16,E16.6)') 'Est=',Est
      write(*,'(A16,E16.6)') 'Egt=',Egt
      write(*,'(A16,E16.6)') 'Elt=',Elt
      write(*,'(A16,E16.6)') 'Eht=',Eht
      write(*,'(A16,E16.6)') 'Est_snow=',Est_snow
      write(*,'(A16,E16.6)') 'Emt_snow=',Emt_snow
      write(*,'(A16,E16.6)') 'Elt_snow=',Elt_snow
      write(*,'(A16,E16.6)') 'Egt_liq =',Egt_liq
      write(*,'(A16,E16.6)') 'Egt_snow=',Egt_snow
      write(*,'(A16,E16.6)') 'Ett=',Ett

      if(forback==1) then
        dTcdt = (Q - Qr)/Qtt
        EJ    = Ett*dTcdt + Er - Ea - Ealpha
        write(*,'(A16, E16.6, I16)') 'time (Q fixed) = ', time, i
      elseif(forback==2) then

        if((sol==0) .or. (sol==1)) then              !Q specified
          dTcdt = (Q - Qr)/Qtt
          EJ    = Ett*dTcdt + Er - Ea - Ealpha
          write(*,'(A16, E16.6, I16)') 'time (Q fixed) = ', time, i

        elseif(sol==2) then                          !OH specified
          dTcdt = (EJ + Ea + Ealpha - Er)/Ett
          Q     = Qtt*dTcdt + Qr
          write(*,'(A16, E16.6, I16)') 'time (E fixed) = ', time, i

        elseif(sol==3) then                          !Q post IC formation; OH otherwise
          if(ri .ne. 0.d0) then
            dTcdt = (Q - Qr)/Qtt
            EJ    = Ett*dTcdt + Er - Ea - Ealpha
            qq = 1
            Efirst=EJ                                !Track EJ through IC formation
            write(*,'(A16, E16.6, I16)') 'time (Q fixed) = ', time, i
          else
            if(qq==1) then                           !No IC for 1st time
              dTcdt = (Q - Qr)/Qtt                   !Use last entropy Efirst
              Efirst  = Ett*dTcdt + Er - Ea - Ealpha
              qq=2                                   !Don't come in here again
            endif
            EJ = Efirst                              !EJ fixed pre IC formation
            dTcdt = (EJ + Ea  + Ealpha - Er)/Ett
            Q     = Qtt*dTcdt + Qr
            write(*,'(A16, E16.6, I16)') 'time (E fixed) = ', time, i
          endif

         endif

        endif

  end subroutine core_en_ent_calc

! If forback = 1, forward; 2=backward
  subroutine get_time(i,time,dt,forback)

    integer          :: i, forback
    double precision :: dt, time

    if(forback==1)  time  = real(i)*dt/1e6
    if(forback==2)  time  = (ttot/1e6) - real(i)*dt/1e6

  end subroutine get_time

  subroutine terminate(ri, rs, flag)

    double precision, intent(in)    :: ri, rs
    integer         , intent(inout) :: flag

    if( (ri .ne. 0.d0) .and. (float(ceiling(rs)) .ne. rc) ) then
      write(*,*) "STOP: IC and Snow layer => Conc incorrect"
      flag=-1
    endif

    if(  (ri .ne. 0.d0) .and. (ri .ge. rs)  ) then 
      write(*,*) "STOP: Inner core occupies whole core!"
      flag=-1
    endif

    if( (float(ceiling(rs)) .ne. rc) .and. (float(ceiling(rs)) .le. ri) ) then
      write(*,*) "STOP: Snow layer occupies whole core!"
      flag=-1
    endif

  end subroutine terminate

  end module core_evolution
