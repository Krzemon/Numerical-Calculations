# Zapisywanie wykresu do pliku
set terminal png

# Ustawienia wykresu
set xlabel "x"
set ylabel "gęstość zliczeń"
set xrange [0:1]
set yrange [0.8:1.15]

# Szerokość binu
binwidth = 1.0 / 10
# Funkcja określająca numer binu
bin(x, width) = width * floor(x / width)
set boxwidth binwidth

set output './plot/markova_a.png'

# Rysowanie histogramu
plot './data/markova_a.dat' using ($1-binwidth/2):2 smooth frequency with boxes title "Łańcuch Markova Δ = 0.5", \
     './data/theory.dat' with lines title "Funkcja gęstości"
