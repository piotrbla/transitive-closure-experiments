printf "26\n30\n14\n" > tile.sizes
./polycc --tile --parallel NusValidation.cpp
gcc -o nus14 -O2 NusValidation.cpp.pluto.c -fopenmp -lm
cp nus14 exp
./exp/nus14
