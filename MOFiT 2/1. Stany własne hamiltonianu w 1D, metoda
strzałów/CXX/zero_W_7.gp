set terminal pngcairo size 1200,800 enhanced font 'Helvetica,16'
set output 'plot.png'
set xlabel 'E [meV]'
set ylabel '$\\psi_{N}$'
set grid
plot 'data.txt' matrix with lines title 'Data 1', 'data.txt' matrix with lines title 'Data 2'
