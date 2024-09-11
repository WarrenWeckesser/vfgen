set -e

echo "Testing scilab"
vfgen scilab ../vf/linearosc.vf
scilab-cli -nwni -noatomsautoload -nb -f test_scilab.sce
