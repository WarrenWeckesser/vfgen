!
! test_dde_solver.f90
!
! Fortran 90 program that will use DDE_SOLVER_M to solve the DDEs defined
! in the vector field file sdd.vf.
!
! If the values computed are not close enough to the known exact solution
! this program will call
!     stop 'FAIL'
! A test script can detect the nontrivial output as a condition for
! failure of the test.
!

program test_dde_solver

use DEFINE_sdd_DDEs
use DDE_SOLVER_M

implicit none

integer, dimension(2) :: NVAR = (/NEQN,NLAGS/)

type(DDE_SOL) :: SOL
type(DDE_OPTS) :: OPTS

double precision, dimension(2) :: TSPAN

integer :: I,J
double precision :: relerr, abserr, stoptime
double precision :: t, y, exact, e
character(100) :: F


! Set the solver parameters: relative error, abs. error, stop time
relerr = 1D-11
abserr = 1D-14
stoptime = 10.0

TSPAN(1) = 0.0
TSPAN(2) = stoptime
OPTS = DDE_SET(RE=relerr,AE=abserr)

SOL = DDE_SOLVER(NVAR,sdd_ddes,sdd_beta,sdd_history,TSPAN,OPTIONS=OPTS)

e = dexp(1.0D0)

do I = 1, SOL%NPTS
    t = SOL%T(I)
    y = SOL%Y(I,1)
    if (t .lt. (e - 1.0D0)) then
        exact = t + 1.0D0
    else if (t .lt. (e*e - 1.0D0)) then
        exact = dexp((t + 1.0D0)/e)
    else
        exact = (e / (3.0D0 - dlog(t + 1.0D0)))**e
    end if
    if (dabs(y - exact) .gt. 1.0D-9) then
        stop 'FAIL'
    end if
end do

end program test_dde_solver
