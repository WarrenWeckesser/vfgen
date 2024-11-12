set -e

echo "Running vfgen taylor:order=4 one_over_x.vf"
vfgen taylor:order=4 one_over_x.vf

echo "Building test_taylorf"
gcc -c one_over_x_taylor4.c
gcc -c test_taylor.c
gcc -o test_taylor test_taylor.o one_over_x_taylor4.o -lm

echo "Running test_taylor"
./test_taylor
