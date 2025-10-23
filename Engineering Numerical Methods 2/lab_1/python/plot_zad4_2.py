#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

if len(sys.argv) != 3:
    print("Uzycie: python3 plot_deltaR_zad4.py <plik_danych> <plik_wykresu>")
    sys.exit(1)

data_file = Path(sys.argv[1])
plot_file = Path(sys.argv[2])

if not data_file.exists():
    print(f"Blad: plik danych '{data_file}' nie istnieje.")
    sys.exit(1)

data = np.loadtxt(data_file, skiprows=1)
r = data[:, 0]
deltaR = data[:, 1:]

plt.figure(figsize=(8, 5))
for i in range(deltaR.shape[1]):
    plt.plot(r, deltaR[:, i], label=f"ΔR_p{i+1}")

plt.xlabel("r")
plt.ylabel("ΔR(r) = R_dok(r) - R_num(r)")
plt.title("Błąd globalny rozwiązania metodą Numerova")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(plot_file, dpi=200)
plt.close()