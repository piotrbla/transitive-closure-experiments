#!/bin/sh

gcc -DZMAXSIZE=900 -o bstexe_900 -O2 bst.cpp -fopenmp -lm
gcc -DZMAXSIZE=1000 -o bstexe_1000 -O2 bst.cpp -fopenmp -lm
gcc -DZMAXSIZE=1100 -o bstexe_1100 -O2 bst.cpp -fopenmp -lm
gcc -DZMAXSIZE=1200 -o bstexe_1200 -O2 bst.cpp -fopenmp -lm
gcc -DZMAXSIZE=1300 -o bstexe_1300 -O2 bst.cpp -fopenmp -lm
gcc -DZMAXSIZE=1400 -o bstexe_1400 -O2 bst.cpp -fopenmp -lm
gcc -DZMAXSIZE=1500 -o bstexe_1500 -O2 bst.cpp -fopenmp -lm
gcc -DZMAXSIZE=1600 -o bstexe_1600 -O2 bst.cpp -fopenmp -lm
gcc -DZMAXSIZE=1700 -o bstexe_1700 -O2 bst.cpp -fopenmp -lm
gcc -DZMAXSIZE=1800 -o bstexe_1800 -O2 bst.cpp -fopenmp -lm
gcc -DZMAXSIZE=1900 -o bstexe_1900 -O2 bst.cpp -fopenmp -lm
gcc -DZMAXSIZE=2000 -o bstexe_2000 -O2 bst.cpp -fopenmp -lm
