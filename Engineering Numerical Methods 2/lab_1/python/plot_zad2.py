#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

if len(sys.argv) != 3:
    print("Uzycie: python3 plot_zad2.py <plik_danych> <plik_wykresu>")
    sys.exit(1)

data_file = Path(sys.argv[1])
plot_file = Path(sys.argv[2])

if not data_file.exists():
    print(f"Blad: plik danych '{data_file}' nie istnieje.")
    sys.exit(1)

data = np.loadtxt(data_file, skiprows=1, delimiter=None)
E_num = data[:,0]
E_theor = data[:,1]
p_points = np.arange(1, len(E_num)+1)  # kolejne punkty p

plt.figure(figsize=(8,5))
plt.plot(p_points, E_theor, "o-", label="E_teor")
plt.plot(p_points, E_num, "s--", label="E_num")
plt.xlabel("Numer punktu p")
plt.ylabel("E")
plt.title("Porownanie energii numerycznej i teoretycznej (l=0)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(plot_file, dpi=200)
plt.close()