#
# test_r.R
#

# options(repos=structure(c(CRAN="https://cloud.r-project.org")))
# install.packages("deSolve")

library(deSolve)

# Load the vector field definition.
source("delayed_logistic.R")


# --- Parameters ---
parameters = c(
    r = 1,
    K = 1,
    tau = 1,
    x0 = 1.0E-4
)

# --- Initial conditions ---
state = c(
    x = 1.0E-4
)

# --- Time values ---
times = c(0, 250)

# --- Call the DDE solver ---
sol = dede(y = state, times = times, func = delayed_logistic, parms = parameters,
           atol = 1e-14, rtol = 1e-12)

# --- Check that the final value is correct.
finalx = sol[2, "x"]
abserr = abs(finalx - 1.0)
error_status = !(abserr < 1e-10)
if (error_status) {
    sprintf("final value is not close to 1; abserr is %13.8e", abserr)
}
quit(status = error_status)
