set -e

echo "Testing julia"
vfgen julia ../vf/linearosc.vf
julia test_julia.jl
