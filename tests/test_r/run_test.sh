
echo "Installing R and the desolve package"
# update-alternatives --install /usr/bin/gfortran gfortran /usr/local/bin/gfortran 999
apt-get -yqq install r-base-core r-cran-desolve

echo "Testing R"
vfgen r ../vf/linearosc.vf

R -f test_r.R --vanilla --slave
