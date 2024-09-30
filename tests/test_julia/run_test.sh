set -e

echo "Generating from linearosc.vf"
vfgen julia ../vf/linearosc.vf

echo "Running test_linearosc.jl"
julia test_linearosc.jl

echo "Generating from linearoscp.vf"
vfgen julia:func=yes ../vf/linearoscp.vf

echo "Running test_linearoscp.jl"
julia test_linearoscp.jl

echo "Generating from sdd.vf"
vfgen julia ../vf/sdd.vf

echo "Running test_sdd.jl"
julia test_sdd.jl
