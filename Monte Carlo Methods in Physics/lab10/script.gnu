# Ustawienie parametrów wykresu
set terminal pngcairo enhanced font 'Verdana,12'
set output 'plots/plot_lab9_t10_n4.png'

# Tytuł wykresu i etykiety osi
set title 'Porównanie wyników numerycznych z dokładnymi dla u(x,t)'
set xlabel 'x'
set ylabel 'u(x,t)'

# Odczyt danych z pliku
plot 'data/lab9_t50_n5.dat' using 1:2 with lines title 'Numeryczne', \
     'data/lab9_t50_n5.dat' using 1:3 with lines title 'Dokładne'
