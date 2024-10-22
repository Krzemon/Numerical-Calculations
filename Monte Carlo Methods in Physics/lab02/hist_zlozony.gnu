# Zapisywanie wykresu do pliku
set terminal png

# Ustawienia wykresu
set xlabel "x"
set ylabel "gęstość zliczeń"
set xrange [0:1]
set yrange [0.:1.15]

# Szerokość binu
binwidth = 1.0 / 10
# Funkcja określająca numer binu
bin(x, width) = width * floor(x / width)
set boxwidth binwidth

set output './plot/zlozony.png'

# Rysowanie histogramu
plot './data/zlozony.dat' using ($1-binwidth/2):2 smooth frequency with boxes title "Rozkład złożony", \
     './data/theory.dat' with lines title "Funkcja gęstości"
