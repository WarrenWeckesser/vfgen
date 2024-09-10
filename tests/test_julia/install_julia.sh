set -e

echo "Installing Julia"
curl -fsSL https://install.julialang.org | sh -s -- -y

echo "Installing Julia DifferentialEquations package"
julia -e 'import Pkg; Pkg.add("DifferentialEquations")'
