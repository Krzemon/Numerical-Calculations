## Rozwiązanie równania Poissona przy użyciu MES 2D z funkcjami kształtu Hermite’a

# Wykonanie:
#### -> cmake -S . -B build
#### -> cmake --build build

### lub 

#### -> mkdir build 
#### -> cd build
#### -> cmake ..
#### -> make

# Program korzysta z zewnętrzych bibliotek:
#### · GSL (GNU Scientific Library)

# Katalogi:
#### · build (należy stworzyć) - Pliki wykonywalne
#### · doc - PDF z poleceniem
#### · src - Pliki źródłowe
#### · include - Pliki nagłówkowe
#### · python - Skrypty do tworzenia wykresów

### ...tworzone później
#### · data - Pliki z wynikami
#### · plots - Pliki z wykresami

· plot_heatmap.py - Skrypt do tworzenia wykresów mapy ciepła macierzy sztywności S
· plot_zad2.py - Skrypt do tworszenia obrazu siatki obliczeniowej elementy, węzły lokalne i globalne
· plot_zad8.py - Skrypt do tworzenia wykresów funkcji rozwiązania u(x,y) analitycznej i numerycznej
· plot_zad8.ipynb - Notatnik Jupyter do tworzenia wykresów funkcji rozwiązania u(x,y) analitycznej i numerycznej z wykorzystaniem biblioteki plotly dla wykresów 3D