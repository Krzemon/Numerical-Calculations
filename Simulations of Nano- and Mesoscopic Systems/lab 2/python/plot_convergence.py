#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
 
# ===================== STYL WYKRESOW =====================
rc("text.latex", preamble=r"\usepackage{lmodern} \usepackage{physics}")
 
plt.rcParams.update({
    "font.size": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "axes.labelsize": 22,
    "axes.linewidth": 1.2,
    "lines.linewidth": 2.5,
    "figure.figsize": (10, 6),
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "text.usetex": True,
    "font.family": "Computer Modern Roman",
})
 
 
# ===================== WCZYTYWANIE DANYCH =====================
def load_data(filename):
    """Wczytuje dwie kolumny: iteracja, E [eV]."""
    data = np.loadtxt(filename, comments='#')
    return data[:, 0], data[:, 1]
 
 
# ===================== RYSOWANIE WYKRESU =====================
def plot_convergence(iterations, energies, output_file):
    fig, ax = plt.subplots()

    # linia ciągła
    ax.plot(iterations, energies, '-', color='steelblue', label='Interpolacja')

    # punkty danych
    ax.plot(iterations, energies, 'o', color='darkred', markersize=4, label='Punkty')

    ax.set_xlabel(r"Iteracja $k$")
    ax.set_ylabel(r"$E^{(k)}\;[\mathrm{eV}]$")
    ax.set_title(r"Zbieznosc metody czasu urojonego")

    ax.grid(True, linestyle="--", alpha=0.3)
    ax.legend()

    fig.savefig(output_file, dpi=300)
    plt.close(fig)

    print(f"Zapisano wykres: {output_file}")
 
 
# ===================== MAIN =====================
def main():
    if len(sys.argv) != 3:
        print("Uzycie: python3 plot_convergence.py <plik_danych> <plik_wyjsciowy>")
        return
 
    data_file   = Path(sys.argv[1])
    output_file = Path(sys.argv[2])
 
    if not data_file.exists():
        print("Plik nie istnieje:", data_file)
        return
 
    try:
        iterations, energies = load_data(data_file)
        plot_convergence(iterations, energies, output_file)
    except Exception as e:
        print(f"Blad: {e}")
        raise
 
 
if __name__ == "__main__":
    main()
    