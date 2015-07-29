#!/usr/bin/gnuplot
set terminal postscript eps font 'Helvetica'
set output
unset border
unset key
set size 0.5,0.5
set zeroaxis lt -1 lw 2
set xtics axis -0.5,0.5,1.5
set ytics axis -1,1,1
quadra(x) = x*x-x-1
set label 1 at 0,-1 point pt 7 ps 1.5
set label 2 at -0.62,0 point pt 6 lw 1.5 ps 1.5
set label 'p' at 2.0,0.3
set label "C(p)U" at 0.1,1.5
set arrow 1 from 2.0,0 to 2.2,0 filled size 0.1,15,60 lw 2
set arrow 2 from 0,1.5 to 0,1.7 filled size 0.1,15,60 lw 2
set arrow 3 from 0,-1 to -0.2,-0.8 filled size 0.1,15,60 lw 1
plot [-1:2] [-1.5:1.7]  quadra(x)

