#
# test_r.R
#

install.packages("deSolve")

library(deSolve)

# Load the vector field definition and the jacobian.
source("linearosc.R")


# --- Initial conditions ---
state = c(x = 1.0, y = 0.0)

# --- Time values ---
times = c(0, pi)

# --- Call the ODE solver ---
sol = ode(y = state, times = times, func = linearosc,
          jactype = "fullusr", jacfunc = linearosc_jac,
          atol = 1e-12, rtol = 1e-8)

# --- Check that the final value is correct.
expected = c(-1, 0)
error = !all(abs(sol[2, c("x", "y")] - expected) < 1e-6)
quit(status = error)
