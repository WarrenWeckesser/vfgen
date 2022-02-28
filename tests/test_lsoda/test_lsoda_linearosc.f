c
c test_lsoda_linearosc.f
c
c Fortran 77 program that uses LSODA to solve the differential equations
c defined in the vector field 'linearosc'.
c
c
      program linearosc

      implicit none

      external linearosc_rhs
      external linearosc_jac

      double precision atol_, rtol_, y_, t_, tout_, tfinal_, rwork_
      integer iwork_
      dimension y_(2), rwork_(52), iwork_(22)
      integer neq_, i_, j_, nsteps_
      integer itol_, iopt_, itask_, istate_, jt_, lrw_, liw_
      double precision x, y
c     --- t range ---
      t_ = 0.0D0
      tfinal_  = 3.1415926535897932384626433D0
      nsteps_ = 100
c     --- Initial conditions ---
      x = 1.0000000000000000D+00
      y = 0.0000000000000000D+00
      y_(1) = x
      y_(2) = y
c     --- Solver tolerances ---
      rtol_ = 1.0D-9
      atol_ = 1.0D-12
      itol_ = 1
c     --- Other LSODA parameters ---
      neq_ = 2
      itask_ = 1
      istate_ = 1
      iopt_ = 0
      lrw_ = 52
      liw_ = 22
c
c     LSODE and LSODA take the same arguments, so either may
c     be used in the loop below.  jt_ must be set as follows:
c     jt_ =  1 for LSODA
c     jt_ = 10 for LSODE, non-stiff (Adams) method
c     jt_ = 21 for LSODE, stiff (BDF) method
c     See the documentation in the Fortran file for more details.
      jt_ = 1

c     --- Call DLSODA in a loop to compute the solution ---
      do 40 i_ = 1,nsteps_
          tout_ = (i_*tfinal_)/nsteps_
          call DLSODA(linearosc_rhs, neq_, y_, t_, tout_,
     &           itol_, rtol_, atol_, itask_, istate_, iopt_,
     &           rwork_, lrw_, iwork_, liw_,
     &           linearosc_jac, jt_)
          if (istate_ .lt. 0) goto 80
40    continue
      x = y_(1)
      y = y_(2)
      if ((abs(x + 1.0D0) .gt. 1D-9) .or. (abs(y) .gt. 1D-9)) then
          stop 'FAIL'
      endif
      stop
80    write (6,89) istate_
89    format(1X,"Error: istate=",I3)
      stop 'FAIL' 
      end
