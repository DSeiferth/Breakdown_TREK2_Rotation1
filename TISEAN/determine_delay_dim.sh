#!/bin/bash
#l=1000000 #number of data to use
b=16 #number of boxes for the partition
D=1000 #maximal time delay
#./mutual ../../SUBSTATES_positive/80.txt -x10 -c2 -b$b -D$D -o mutual_+80.dat

#./mutual ../../SUBSTATES_positive/80.txt -x10 -c2 -b200 -D$D -o mutual_+80b200.dat

gnuplot -persist <<-EOFMarker
    #set title sprintf('Embedding: m=%d tau=%d', $m, $d) font ",14" 
    set xlabel "delay tau" font ",14"
    set ylabel "mutual information" font ",14"
    set key right bottom font ",14"
    set logscale xy
    
    plot 'mutual_+80.dat' with lines title "+80 mV", 'mutual_+80b200.dat' with lines title "+80 mV b200"
    #set terminal postscript enhanced color solid
    #set output "mutual_information.eps"
    #set encoding iso_8859_15
    #replot
    #set terminal x11
EOFMarker
