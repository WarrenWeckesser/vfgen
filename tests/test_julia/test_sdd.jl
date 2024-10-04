using DifferentialEquations
using Printf

include("sdd.jl")

function history(p, t; idxs = nothing)
    if typeof(idxs) <: Number
        1.0
    else
        [1.0]
    end
end

function dependent_lag(u, p, t)
    t - (log(u[1]) - 1.0)
end

u0 = [1.0]
tfinal = 10.0
tspan = (0.0, tfinal)

prob = DDEProblem(sdd!, u0, history, tspan;
                  dependent_lags = [dependent_lag])
alg = MethodOfSteps(Tsit5())
sol = solve(prob, alg, abstol=5e-15, reltol=1e-13)
xfinal = sol[end][1]

e = Base.MathConstants.e
exact = (e / (3.0 - log(tfinal + 1.0)))^e

relerr = abs(xfinal - exact)/exact
if relerr > 2e-10
    print("Final value ", xfinal, " is not close to the exact solution ", exact)
    @printf("; rel err is %.4e\n", relerr)
    exit(-1)
end
