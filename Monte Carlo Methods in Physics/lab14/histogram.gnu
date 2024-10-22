set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'plots/histogram.png'
set title "Histogram of Sampled Points"
set xlabel "r"
set ylabel "Frequency"
set style data histograms
set style fill solid 1.0 border -1
plot 'data/histogram.dat' using 2:xtic(1) title 'Sampled Points'
