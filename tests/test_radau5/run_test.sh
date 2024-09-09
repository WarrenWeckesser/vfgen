
vfgen radau5 ../vf/linearosc.vf

gfortran -c radau5.f
gfortran -c dc_lapack.f
gfortran -c linearosc_rhs.f
gfortran -c test_linearosc.f
gfortran test_linearosc.o linearosc_rhs.o radau5.o dc_lapack.o -llapack -o test_linearosc

./test_linearosc
