#!/usr/bin/env gnuplot

set view equal xyz
set xlabel "x"
set ylabel "y"
set zlabel "z"
set xrange [-0.2:2]
splot "a.dat" i 0 t "corner", \
      "a.dat" i 1 t "collocation point", \
      "a.dat" i 2 usi 1:2:3:4:5:6 w vec t "normal", \
      "a.dat" i 3 usi 1:2:3:4:5:6 w vec t "v at collocation point", \
      "a.dat" i 7 usi 1:2:3:($4*5):($5*5):($6*5) w vec t "induced v x10" lc rgb "red"
