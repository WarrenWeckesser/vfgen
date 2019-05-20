#
# test_octave.m
#

# Load the vector field definition and the jacobian.
source "linearosc.m";


# --- Initial conditions ---
x(1) = 1.0;
x(2) = 0.0;

# --- Time values ---
t0 = 0.0;
t1 = pi;
numpoints = 201;
t = linspace(t0, t1, numpoints);

# --- Solver error tolerances ---
lsode_options("relative tolerance", 1e-8);
lsode_options("absolute tolerance", 1e-12);

# --- Call the ODE solver ---
xsol = lsode({"linearosc_vf", "linearosc_jac"}, x, t);

final = xsol(end, :);
expected = [-1 0];
pass = max(abs(final - expected)) < 1e-6;

exit(~pass)
