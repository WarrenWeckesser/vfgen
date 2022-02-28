c
c test_lsoda_pidecay.f
c
c Fortran 77 program that uses LSODA to solve the differential equation
c defined in the vector field 'pidecay'.  The differential equation is
c
c   x' = -Pi*x
c
c The exact solution is
c
c   x(t) = x(0)*exp(-Pi*t)
c
      program linearosc

      implicit none

      external pidecay_rhs
      external pidecay_jac

      double precision atol, rtol, y, t, tout, tfinal, rwork
      integer iwork
      dimension y(1), rwork(52), iwork(22)
      integer neq, i, j, nsteps
      integer itol, iopt, itask, istate, jt, lrw, liw

c     --- t range ---
      t = 0.0D0
      tfinal  = 2.0D0
      nsteps = 5

c     --- Initial conditions ---
      y(1) = 1.0

c     --- Solver tolerances ---
      rtol = 1.0D-9
      atol = 1.0D-12
      itol = 1

c     --- Other LSODA parameters ---
      neq = 1
      itask = 1
      istate = 1
      iopt = 0
      lrw = 52
      liw = 22

c     LSODE and LSODA take the same arguments, so either may
c     be used in the loop below.  jt must be set as follows:
c     jt =  1 for LSODA
c     jt = 10 for LSODE, non-stiff (Adams) method
c     jt = 21 for LSODE, stiff (BDF) method
c     See the documentation in the Fortran file for more details.
      jt = 1

c     --- Call DLSODA in a loop to compute the solution ---
      do 40 i = 1, nsteps
          tout = (i*tfinal) / nsteps
          call DLSODA(pidecay_rhs, neq, y, t, tout,
     &           itol, rtol, atol, itask, istate, iopt,
     &           rwork, lrw, iwork, liw,
     &           pidecay_jac, jt)
          if (istate .lt. 0) goto 80
40    continue

c The expected value is
c
c     exp(-pi*tfinal) = exp(-pi*2) = 0.0018674427317079889
      if (abs(y(1) - 0.0018674427317079889D0) .gt. 1D-9) then
          stop 'FAIL'
      endif
      stop
80    write (6,89) istate
89    format(1X,"Error: istate=",I3)
      stop 'FAIL' 
      end
