## Metoda shootingu w 1D, metoda różnic skończonych, metoda Numerowa.

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
------  
---
------
### Deklaracje funkcji z Numerical Recipes, w kodzie cpp tworzymy:
#### · float bessj0(float r);
#### · float bessj1(float r);

### Kompilacja z numerical_recipes:
#### gcc test_bessel.c /opt/NR/numerical_recipes.c -lm -o test_bessel
#### Parametr -lm dodaje bibliotekę matematyczną (libm), bo funkcje NR mogą korzystać np. z sin() czy sqrt().
