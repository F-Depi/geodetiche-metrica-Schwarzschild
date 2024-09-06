# ./main.x 1.8 -0.042 -s 1 -t 1000 -r 2 -B 500 -h 2e-3 -a 0 -f keep/infall2_stability/l1.8E-0.042_h2.0e-03
# ./main.x 1.8 -0.042 -s 1 -t 1000 -r 2 -B 500 -h 1e-3 -a 0 -f keep/infall2_stability/l1.8E-0.042_h1.0e-03

#!/bin/bash

# Common parameters
L=1.8
E=-0.042
S=1
T=1000
R=2
B=500
A=0
DIR="keep/infall2_stability"

# Run the commands with varying 'h' values
for H in 2.0e-03 1.0e-03 8.0e-04 4.0e-04 2.0e-04 1.0e-04 5.0e-05 2.5e-05 1.0e-05; do
    ./main.x $L $E -s $S -t $T -r $R -B $B -h $H -a $A -f "$DIR/l${L}_E${E}_h${H}"
done

