# Ścieżka do pliku danych
data_file = "build/data/temperature.txt"

# Ścieżka do pliku wyjściowego (outputu)
output_file = 'temperature.png'

set datafile separator ","

# Ustawienia wykresu
set title "Temperatura molekuł w rurze w zależności od numeru iteracji"
set xlabel "Numer iteracji"
set ylabel "temperatura"

# Styl punktów z mniejszymi punktami
set style line 1 lt 1 lc rgb "blue" lw 2 pt 7 ps 0.4

# Zapisz wykres do pliku wyjściowego w formacie PNG
set terminal pngcairo enhanced font 'Verdana,12'
set output output_file
plot data_file using 1:2 with points ls 1 notitle
set output
