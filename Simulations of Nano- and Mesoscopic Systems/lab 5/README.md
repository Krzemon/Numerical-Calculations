# Transport elektronowy w układzie 2D: kwantowy kontakt punktowy (QPC)

> **Sprawozdanie z laboratorium znajduje się w katalogu `doc/`**

---

## Opis projektu

Symulacja transportu elektronowego w kwantowym kontakcie punktowym (QPC) metodą QTBM (Quantum Transmitting Boundary Method).

**Realizowane zadania:**
1. Relacja dyspersji E(kx) dla jednorodnego kanału
2. Mody propagujące się dla różnych energii
3. Test poprawności: T + R = liczba modów
4. Przewodność G(Vgates) i mapa potencjału V(x,y)
5. Mapy gęstości |ψ|² dla różnych napięć bramek

---

## Kompilacja i uruchomienie

### Metoda 1:
```bash
cd build && make -j4
```

### Metoda 2:
```bash
cmake -S . -B build
cmake --build build
```

### Metoda 3:
```bash
mkdir -p build
cd build
cmake ..
make
```

---

## Wyniki

- **Dane:** `data/*.dat`
- **Wykresy:** `plots/zad_*/`

---

## Wymagania

- CMake ≥ 3.10
- C++17
- Eigen3
- Python 3 + matplotlib + numpy (do wykresów)
