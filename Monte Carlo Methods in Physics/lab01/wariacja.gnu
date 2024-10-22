set terminal pngcairo
set output "plot/wariancja.png"

set xlabel "Numer iteracji"
set ylabel "Względny błąd wariancji"

# Ustawienie znaczników na osi X w postaci wykładniczej
set format x "10^{%L}"
set format y "10^{%L}"

# Ustawienie zakresu osi Y
set yrange [0.0000001:1]

# Ustawienie podziałki osi X i Y na logarytmiczną
set logscale y 10
set logscale x 10

# Rysowanie trzech krzywych z różnych plików
plot "./data/p01.dat" using (10**$1):3 with linespoints linecolor rgb "blue" pointtype 7 pointsize 1.5 title "p=0.1", \
     "./data/p05.dat" using (10**$1):3 with linespoints linecolor rgb "red" pointtype 7 pointsize 1.5 title "p=0.5", \
     "./data/p09.dat" using (10**$1):3 with linespoints linecolor rgb "green" pointtype 7 pointsize 1.5 title "p=0.9"
