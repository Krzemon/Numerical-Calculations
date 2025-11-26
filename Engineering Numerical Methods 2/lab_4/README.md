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


## Zadanie 1: 
· plot_zad1.py – mapa ciepła macierzy S.

## Zadanie 2:
· plot_zad2_E.py – wykres wartości własnych E_mu(alpha) dla podanego pliku danych, obsługuje dowolne M.
· plot_zad2_psi.py – wykres funkcji własnych u_mu(x) dla α=1.4 dla dowolnego M.

## Zadanie 3: 
· plot_zad3.py – wykres wartości własnych dla M=10 i M=30, obsługujący bloki danych w jednym pliku.