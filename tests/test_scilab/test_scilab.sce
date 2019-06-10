//
// test_scilab.sce
//

// Load the vector field definition and the jacobian.
exec('linearosc.sci');


// --- Initial conditions ---
x0 = [1.0; 0.0];

// --- Time values ---
t0 = 0.0;
t1 = %pi;
numpoints = 201;
t = linspace(t0, t1, numpoints);

// --- Call the ODE solver ---
sol = ode(x0, t0, t, 1e-9, 1e-12, linearosc_vf, linearosc_jac);


final = sol(:, $);
expected = [-1; 0];
pass = max(abs(final - expected)) < 1e-6;
exit(pass - 1);
