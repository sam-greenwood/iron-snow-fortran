  program solid_fraction

    implicit none

    double precision, parameter :: pi = 3.1415926535897931d0
    integer         , parameter :: nt=1000, ns=1000
    double precision   :: cl, cl0, V, Vc, phi_t
    double precision   :: r(ns), phi(ns), int1(ns), int2(ns)
    integer :: i, j

    open(10,file="fort.888")
    open(11,file="phi_ts")
   
    
    do i = 1, nt
      phi_t = 0.d0; int1=0.d0
      do j = 1, 51
        read(10,*,end=100) r(j), cl, cl0, phi(j) !phi=0 here
        print *, phi(j)
      enddo
      do j = 1, ns
        read(10,*,end=100) r(j), cl, cl0, phi(j)
        if(phi(j)==0.d0) stop "phi should not be zero!"

        int1(j) = r(j)*r(j)*phi(j)
        int2(j) = r(j)*r(j)
      enddo
      read(10,*)
      read(10,*)

      Vc = 4.d0 * pi * (r(ns)**3) / 3.d0
      print *, r(ns), Vc

      call splat(int1,r,ns)
      phi_t = 4.d0 * pi * int1(ns-1)  
      call splat(int2,r,ns)
      V     = 4.d0 * pi * int2(ns-1)

      write(11,'(I16,4E16.8)') i, V, phi_t, phi_t/Vc
    enddo

100 continue

    close(10)
    close(11)

  end program solid_fraction
