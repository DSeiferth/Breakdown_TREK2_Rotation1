#!/bin/bash
l=100000 #number of data to use
m=1 #minimal embedding dimensions of the vectors
# -M#,# # of components, max. embedding dimension of the vectors
# default 1,5 
M=10 
d=5 #	delay of the vectors 
f=10 #	ratio factor (default 2)

# This program looks for the nearest neighbors of all data points in m dimensions and iterates these neighbors one step (more precisely delay steps) into the future.
# f the ratio of the distance of the iteration and that of the nearest neighbor exceeds a given threshold the point is marked as a wrong neighbor. 
# We implemented a new second criterion. If the distance to the nearest neighbor becomes smaller than the standard deviation of the data devided by the threshold, the point is omitted.
# This turns out to be a stricter criterion, but can show the effect that for increasing embedding dimensions the number of points which enter the statistics is so small, that the whole statistics is meanlingless.

./false_nearest ../../SUBSTATES/160.txt -l$l -x10 -c2 -M$m,$M -d$d -f$f -o ffn-160_tau5.dat
./false_nearest ../../SUBSTATES/140.txt -l$l -x10 -c2 -M$m,$M -d$d -f$f -o ffn-140_tau5.dat
./false_nearest ../../SUBSTATES/120.txt -l$l -x10 -c2 -M$m,$M -d$d -f$f -o ffn-120_tau5.dat
./false_nearest ../../SUBSTATES/100.txt -l$l -x10 -c2 -M$m,$M -d$d -f$f -o ffn-100_tau5.dat
./false_nearest ../../SUBSTATES/80.txt -l$l -x10 -c2 -M$m,$M -d$d -f$f -o ffn-80_tau5.dat
./false_nearest ../../SUBSTATES/60.txt -l$l -x10 -c2 -M$m,$M -d$d -f$f -o ffn-60_tau5.dat
./false_nearest ../../SUBSTATES/40.txt -l$l -x10 -c2 -M$m,$M -d$d -f$f -o ffn-40_tau5.dat
./false_nearest ../../SUBSTATES/20.txt -l$l -x10 -c2 -M$m,$M -d$d -f$f -o ffn-20_tau5.dat

gnuplot -persist <<-EOFMarker
    set title sprintf('100k datapoints, embedding delay: tau=%d', $d) font ",14" 
    set xlabel "embedding dimension m" font ",14"
    set ylabel "fraction of false nearest neighbours" font ",14"
    set key right top font ",14"
    #set xrange [0: 0.1]
    #set xrange [3: 10]
    #set logscale xy
    
    plot 'ffn-160_tau5.dat' with lines title "-160 mV", 'ffn-140_tau5.dat' with lines title "-140 mV", 'ffn-120_tau5.dat' with lines title "-120 mV", 'ffn-100_tau5.dat' with lines title "-100 mV", 'ffn-80_tau5.dat' with lines title "-80 mV", 'ffn-60_tau5.dat' with lines title "-60 mV", 'ffn-40_tau5.dat' with lines title "-40 mV", 'ffn-20_tau5.dat' with lines title "-20 mV"
    set terminal postscript enhanced color solid
    set output "FalseNearestNeighbour_tau=5.eps"
    set encoding iso_8859_15
    replot
    set terminal x11
EOFMarker

gnuplot -persist <<-EOFMarker
    set title sprintf('100k datapoints, embedding delay: tau=%d', $d) font ",14" 
    set xlabel "embedding dimension m" font ",14"
    set ylabel "fraction of false nearest neighbours" font ",14"
    set key right top font ",14"
    set xrange [0: 0.1]
    set xrange [3: 10]
    #set logscale xy
    
    plot 'ffn-160_tau5.dat' with lines title "-160 mV", 'ffn-140_tau5.dat' with lines title "-140 mV", 'ffn-120_tau5.dat' with lines title "-120 mV", 'ffn-100_tau5.dat' with lines title "-100 mV", 'ffn-80_tau5.dat' with lines title "-80 mV", 'ffn-60_tau5.dat' with lines title "-60 mV", 'ffn-40_tau5.dat' with lines title "-40 mV", 'ffn-20_tau5.dat' with lines title "-20 mV"
    set terminal postscript enhanced color solid
    set output "FalseNearestNeighbour_tau=5_zoom.eps"
    set encoding iso_8859_15
    replot
    set terminal x11
EOFMarker

