#
# test_matlab_with_octave.m
#

# Initial conditions
x(1) = 1.0;
x(2) = 0.0;

# Time values
t0 = 0.0;
t1 = pi;
numpoints = 201;
t = linspace(t0, t1, numpoints);

# Solver settings
opts = odeset("AbsTol", 1e-12, "RelTol", 1e-8, "Jacobian", @linearosc_jac);

# Call the ODE solver
[t, xsol] = ode45(@linearosc_vf, t, x, opts);

final = xsol(end, :);
expected = [-1 0];
pass = max(abs(final - expected)) < 1e-6;

exit(~pass)
