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
def load_map(filename):
    """
    Wczytuje plik z blokami (jeden blok = jeden F).
    Zwraca siatke: F_vals, t_vals, Z (macierz p1).
    """
    data = np.loadtxt(filename, comments='#')
    # kolumny: F, t, p1
    F_all = data[:, 0]
    t_all = data[:, 1]
    p1_all = data[:, 2]
 
    F_vals = np.unique(F_all)
    t_vals = np.unique(t_all)
 
    Z = np.zeros((len(F_vals), len(t_vals)))
 
    F_idx = {v: i for i, v in enumerate(F_vals)}
    t_idx = {v: i for i, v in enumerate(t_vals)}
 
    for fi, ti, pi in zip(F_all, t_all, p1_all):
        Z[F_idx[fi], t_idx[ti]] = pi
 
    return F_vals, t_vals, Z
 
 
# ===================== RYSOWANIE =====================
def plot_map(F_vals, t_vals, Z, output_file):
    T, F = np.meshgrid(t_vals, F_vals)
 
    fig, ax = plt.subplots()
    im = ax.pcolormesh(T, F, Z, shading='auto', cmap='viridis', vmin=0.0, vmax=1.0)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"$|\langle\Psi(t)|1\rangle|^2$")
    ax.set_xlabel(r"$t\;[\mathrm{ns}]$")
    ax.set_ylabel(r"$F\;[\mathrm{kV/cm}]$")
    ax.set_title(
        r"Prawdopodobieństwo stanu wzbudzonego "
        r"$|\langle\Psi(t)|1\rangle|^2$"
    )
    fig.savefig(output_file)
    plt.close(fig)
    print(f"Zapisano: {output_file}")
 
 
# ===================== MAIN =====================
def main():
    if len(sys.argv) != 3:
        print("Użycie: python3 plot_evolution_map.py "
                        "<map_p1_t_F.dat> <output.png>")
        return
 
    file_in = Path(sys.argv[1])
    out     = Path(sys.argv[2])
 
    if not file_in.exists():
        print(f"Plik nie istnieje: {file_in}")
        return
 
    try:
        F_vals, t_vals, Z = load_map(file_in)
        plot_map(F_vals, t_vals, Z, out)
    except Exception as e:
        print(f"Błąd: {e}")
        raise
 
 
if __name__ == "__main__":
    main()