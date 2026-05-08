## Wyznaczanie charakterystyki prądowo-napięciowej diody rezonansowo-tunelowej (RTD) oraz zastosowanie przybliżenia adiabatycznego do wyznaczenia zjawiska kwantyzacji konduktancji w kwantowym kontakcie punktowym (QPC).

---

## Sprawozdanie z laboratorium znajduje się w katalogu doc/

# Wykonanie:
#### -> cmake -S . -B build
#### -> cmake --build build

### lub 

#### -> mkdir build 
#### -> cd build
#### -> cmake ..
#### -> make

## Wymagania

- CMake 3.10 lub nowszy
- Kompilator C++ z obsługą C++17
- Biblioteka Eigen3
- OpenMP (opcjonalnie, dla przyspieszenia obliczeń)
- Python 3 z bibliotekami: numpy, matplotlib

### Zadania zrównoleglone (OpenMP)

- **Zadanie 3b (RTD I-V)**: Pętla po napięciu bias
- **Zadanie 4b (QPC G(E))**: Pętla po energii
- **Zadanie 4c (QPC G(Vg))**: Pętla po napięciu bramki

Program domyślnie wykona się na **jednym wątku**, chyba że ustawisz zmienną środowiskową `OMP_NUM_THREADS`.

### Uruchomienie z wieloma wątkami (OpenMP)

Aby wykorzystać wszystkie dostępne rdzenie procesora:

```bash
export OMP_NUM_THREADS=$(nproc)
./nano_mezo_4
```

Lub dla konkretnej liczby wątków (np. 4):

```bash
export OMP_NUM_THREADS=4
./nano_mezo_4
```

Lub w jednej linii:

```bash
OMP_NUM_THREADS=$(nproc) ./nano_mezo_4
```

