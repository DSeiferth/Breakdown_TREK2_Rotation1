#!/bin/bash
l=500000
m=4
d=5
s=50
#./lyap_r ../SUBSTATES/160.txt -l$l -x10 -c2 -m$m -d$d -s$s -o test_160.dat
#./lyap_r ../SUBSTATES/140.txt -l$l -x10 -c2 -m$m -d$d -s$s -o test_140.dat
#./lyap_r ../SUBSTATES/120.txt -l$l -x10 -c2 -m$m -d$d -s$s -o test_120.dat
#./lyap_r ../SUBSTATES/100.txt -l$l -x10 -c2 -m$m -d$d -s$s -o test_100.dat
#./lyap_r ../SUBSTATES/80.txt -l$l -x10 -c2 -m$m -d$d -s$s -o test_80.dat
#./lyap_r ../SUBSTATES/60.txt -l$l -x10 -c2 -m$m -d$d -s$s -o test_60.dat
#./lyap_r ../SUBSTATES/40.txt -l$l -x10 -c2 -m$m -d$d -s$s -o test_40.dat
#./lyap_r ../SUBSTATES/20.txt -l$l -x10 -c2 -m$m -d$d -s$s -o test_20.dat

#./lyap_k ../SUBSTATES/20.txt -l$l -x10 -c2 -m3 -M5 -d$d -s$s -o test_k20.dat

gnuplot -persist <<-EOFMarker
    set title sprintf('Embedding: m=%d, tau=%d, 500k data points', $m, $d) font ",14" 
    set xlabel "t_{evolve}" font ",14"
    set ylabel "S(t_{evolve})" font ",14"
    set key right bottom font ",14"
    f(x) = y0 + m*x
    y0 = -4
    m = 0.01
    #fit f(x) 'test_120.dat' via y0, m
    title_f(y0,m) = sprintf('fit f(t_{evolve}) =%.2f + %.4ft_{evolve} ',  y0, m)
    # plot  f(x) t title_f(y0,m),
    
    plot 'test_160.dat' with lines title "-160 mV", 'test_140.dat' with lines title "-140 mV", 'test_120.dat' with lines title "-120 mV", 'test_100.dat' with lines title "-100 mV", 'test_80.dat' with lines title "-80 mV", 'test_60.dat' with lines title "-60 mV", 'test_40.dat' with lines title "-40 mV", 'test_20.dat' with lines title "-20 mV"
    #, 'test_k20.dat' with lines title "-20 mV Kantz"
    set terminal postscript enhanced color solid
    set output "all_mV_500k_datapoints.eps"
    set encoding iso_8859_15
    replot
    set terminal x11
EOFMarker
