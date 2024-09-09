c
c test_linearosc.f
c
c Fortran 77 program that uses RADAU5 to solve the differential equations
c defined in the vector field 'linearosc'
c
      program linearosc
      implicit none
      integer nd, lwork, liwork
      parameter (nd=2, lwork=60, liwork=26)
      double precision y_, rpar, work, iwork
      dimension y_(nd), work(lwork), iwork(liwork)
      integer n, ipar, ijac, mljac, mujac, imas, mlmas, mumas
      integer iout, itol, idid
      dimension rpar(0)
      integer i
      double precision t, tstop
      double precision atol, rtol, h
      double precision x, y
      external linearosc_rhs
      external linearosc_jac
      external linearosc_out

      n = 2
      ijac = 1
      mljac = n
      imas = 0
      iout = 0

c     --- t range ---
      t = 0.0D0
      tstop  = 3.1415926535897932384626433D0

c     --- Initial conditions ---
      x = 1.0000000000000000D+00
      y = 0.0000000000000000D+00
      y_(1) = x
      y_(2) = y

c     --- Solver tolerances ---
      rtol = 1.0D-14
      atol = 5.0D-15
      itol = 0

c     --- Initial step size ---
      h = 1.0D-5

c     --- Set default values ---
      do i = 1, 20
          iwork(i) = 0
          work(i) = 0.0D0
      end do

c     --- Call RADAU5 ---
      call radau5(n, linearosc_rhs, t, y_, tstop, h,
     &           rtol, atol, itol,
     &           linearosc_jac, ijac, mljac, mujac,
     &           linearosc_rhs, imas, mlmas, mumas,
     &           linearosc_out, iout,
     &           work, lwork, iwork, liwork, rpar, ipar, idid)
      x = y_(1)
      y = y_(2)
      if ((abs(x + 1.0D0) .gt. 1D-10) .or. (abs(y) .gt. 1D-10)) then
          stop 1
      endif
      stop
      end
