
echo "Installing scilab"
apt-get -yqq install scilab-cli
scilab-cli -version

echo "Testing scilab"
vfgen scilab ../vf/linearosc.vf
scilab-cli -nwni -noatomsautoload -nb -f test_scilab.sce -quit
