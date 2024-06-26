#!/bin/bash
#compile program and time
echo "compiling CoOMBE"
make
echo "compile finished"
echo "run CoOMBE (timed)"
time ./CoOMBE < prop_k.dat
#option to remove binary after program run
rm ./CoOMBE
echo "run complete"
