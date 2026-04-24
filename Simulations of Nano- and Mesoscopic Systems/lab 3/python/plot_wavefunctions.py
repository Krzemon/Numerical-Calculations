#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
 
# ===================== STYL =====================
 
plt.rcParams["text.usetex"] = True
plt.rcParams["text.latex.preamble"] = r"""
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage[polish]{babel}
\usepackage{lmodern}
\usepackage{physics}
"""

plt.rcParams.update({
    "font.size": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "axes.labelsize": 22,
    "axes.linewidth": 1.2,
    "lines.linewidth": 2.5,
    "figure.figsize": (10, 7),
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "text.usetex": True,
    "font.family": "Computer Modern Roman",
})
 
# ===================== WCZYTYWANIE =====================
def load_data(filename):
    """Kolumny: x[nm], psi0, psi1."""
    data = np.loadtxt(filename, comments='#')
    return data[:, 0], data[:, 1], data[:, 2]
 
 
# ===================== RYSOWANIE =====================
def plot_wavefunctions(x, psi0, psi1, output_file):
    fig, ax = plt.subplots()
 
    ax.plot(x, psi0, '-',  color='steelblue', label=r"$\psi_0$ (stan wiążący)")
    ax.plot(x, psi1, '--', color='tomato', label=r"$\psi_1$ (stan antywiążący)")
 
    ax.axhline(0, color='black', linewidth=0.8, linestyle=':')
 
    ax.set_xlabel(r"$x\;[\mathrm{nm}]$")
    ax.set_ylabel(r"$\psi(x)\;[\mathrm{a.u.}^{-1/2}]$")
    ax.set_title(r"Funkcje falowe dwóch najniższych stanów ($F=0$)")
    ax.legend(fontsize=14)
    ax.grid(True, linestyle="--", alpha=0.3)
 
    fig.savefig(output_file)
    plt.close(fig)
    print(f"Zapisano: {output_file}")
 
 
# ===================== MAIN =====================
def main():
    if len(sys.argv) != 3:
        print("Użycie: python3 plot_wavefunctions.py "
                        "<wavefunctions.dat> <output.png>")
        return
 
    file_in = Path(sys.argv[1])
    file_out = Path(sys.argv[2])
 
    if not file_in.exists():
        print(f"Plik nie istnieje: {file_in}")
        return
 
    try:
        x, psi0, psi1 = load_data(file_in)
        plot_wavefunctions(x, psi0, psi1, file_out)
    except Exception as e:
        print(f"Błąd: {e}")
        raise
 
 
if __name__ == "__main__":
    main()