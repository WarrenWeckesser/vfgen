using DifferentialEquations

include("linearosc.jl")

u0 = [1.0, 0.0]
tspan = (0.0, pi)

prob = ODEProblem(linearosc!, u0, tspan)
sol = solve(prob, abstol = 1e-14, reltol = 1e-12)
x, y = sol[end]

if (abs(x + 1) > 1e-10) || abs(y) > 1e-10
    println("Final value is not close to (-1, 0).")
    exit(-1)
end
