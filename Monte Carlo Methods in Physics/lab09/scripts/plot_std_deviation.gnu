set term png
set output "../plots/std_deviation.png"
set title "Mapa kolorów odchylenia standardowego potencjału z metody Monte Carlo"
set xlabel "X"
set ylabel "Y"
set palette defined (0 "blue", 1 "red")
set xrange [0:30]
set yrange [0:30]
plot "../data/std_deviation.txt" matrix with image notitle