#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

if len(sys.argv) != 3:
    print("Uzycie: python3 plot_zad3.py <plik_danych> <plik_wykresu>")
    sys.exit(1)

data_file = Path(sys.argv[1])
plot_file = Path(sys.argv[2])

if not data_file.exists():
    print(f"Blad: plik danych '{data_file}' nie istnieje.")
    sys.exit(1)

# Okreslamy l na podstawie nazwy pliku: zad_3_0.dat -> l=0, zad_3_1.dat -> l=1
l = 0 if "zad_3_0" in data_file.stem else 1

data = np.loadtxt(data_file, skiprows=1)
r = data[:,0]
deltaR = data[:,1:]  # kolumny: ΔR_p1, ΔR_p2, ...

plt.figure(figsize=(8,5))
for i in range(deltaR.shape[1]):
    plt.plot(r, deltaR[:,i], label=f"ΔR_p{i+1}")

plt.xlabel("r")
plt.ylabel("ΔR(r) = R_teor(r) - R_num(r)")
plt.title(f"Blad globalny rozwiazania dla czterech pierwszych zer Bessela (l={l})")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(plot_file, dpi=200)
plt.close()
print(f"Zapisano wykres do: {plot_file.resolve()}")