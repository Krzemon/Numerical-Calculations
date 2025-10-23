#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def main():
    if len(sys.argv) != 3:
        print("Uzycie: python3 plot_zad1.py <plik_danych> <plik_wykresu>")
        sys.exit(1)

    data_file = Path(sys.argv[1])
    plot_file = Path(sys.argv[2])

    if not data_file.exists():
        print(f"Blad: plik danych '{data_file}' nie istnieje.")
        sys.exit(1)

    # Wczytaj dane (pomijamy pierwszą linię nagłówka)
    try:
        data = np.loadtxt(data_file, skiprows=1)
    except Exception as e:
        print(f"Błąd przy wczytywaniu danych: {e}")
        sys.exit(1)

    if data.ndim != 2 or data.shape[1] < 2:
        print("Blad: plik danych musi miec co najmniej 2 kolumny.")
        sys.exit(1)

    E, Rn = data[:, 0], data[:, 1]

    plt.figure(figsize=(8, 5))
    plt.plot(E, Rn, label="R_n(E)")
    plt.xlabel("E")
    plt.ylabel("R_n (r = L)")
    plt.title("R_n = f(E) dla l = 0 (R0=R1=1.0)")
    plt.grid(True)
    plt.legend()

    L = 1.0
    zeros_bessel_l0 = [2.4048, 5.5200, 8.6537, 11.7915]
    energies = [0.5 * (a / L) ** 2 for a in zeros_bessel_l0]
    for E_theor in energies:
        plt.axvline(E_theor, color='r', linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig(plot_file, dpi=200)
    plt.close()

    print(f"Zapisano wykres do: {plot_file.resolve()}")

if __name__ == "__main__":
    main()