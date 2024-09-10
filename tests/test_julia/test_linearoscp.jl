using DifferentialEquations

include("linearoscp.jl")

u0 = [1.0, 0.0]
tspan = (0.0, pi)
params = [0.5]

prob = ODEProblem(linearoscp!, u0, tspan, params)
sol = solve(prob, abstol = 1e-14, reltol = 1e-12)
x, y = sol[end]

if (abs(x) > 1e-10) || abs(y - 1) > 1e-10
    println("Final value is not close to (0, 1).")
    exit(-1)
end

sqmag = linearoscp_func(sol[end], params, sol.t[end])
if abs(sqmag - 1) > 1e-10
    println("Final squared magnitude is not close to 1.")
    exit(-1)
end

theta = linearoscp_theta(sol[end], params, sol.t[end])
if abs(theta - 0.5*tspan[end]) > 1e-10
    println("Final theta does not match expected value.")
    println("theta    = ", theta)
    println("expected = ", 0.5*tspan[end])
    exit(-1)
end
