
echo "Installing R and the desolve package"
apt-get -yqq install r-base-core r-cran-desolve

echo "Testing R"
vfgen r ../vf/linearosc.vf

R -f test_r.R --vanilla --slave
