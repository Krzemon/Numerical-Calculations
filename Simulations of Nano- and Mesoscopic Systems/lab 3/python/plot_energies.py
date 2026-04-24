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
    return np.loadtxt(filename, comments='#')
 
 
# ===================== RYSOWANIE =====================
def plot_n_states(data, output_file, n):
    F   = data[:, 0]
    E_table = [data[:, i] for i in range(1, n+1)]
 
    fig, ax = plt.subplots()

    for i in range(n):
        label = r"$E_{}$".format(i)
        if i == 0:
            label += " (wiążący)"
        elif i == 1:
            label += " (antywiążący)"
        color = ['steelblue', 'tomato', 'seagreen', 'darkorange'][i % 4]
        linestyle = '-' if i % 2 == 0 else '--'
        ax.plot(F, E_table[i], linestyle, color=color, label=label)
 
    ax.set_xlabel(r"$F\;[\mathrm{kV/cm}]$")
    ax.set_ylabel(r"$E\;[\mathrm{meV}]$")
    ax.set_title(r"Energie {n} najniższych stanów z polem elektrycznym")
    ax.legend(fontsize=14)
    ax.grid(True, linestyle="--", alpha=0.3)
    fig.savefig(output_file)
    plt.close(fig)
    print(f"Zapisano: {output_file}")
 
# ===================== MAIN =====================
def main():
    if len(sys.argv) != 4:
        print("Użycie: python3 plot_energies.py " 
                        "<states.dat> <out_4.png> <out_2.png>")
        return
 
    file_4 = Path(sys.argv[1])
    out_4 = Path(sys.argv[2]) # wykres 4 stanów
    out_2 = Path(sys.argv[3]) # wykres 2 stanów
 
    if not file_4.exists():
        print(f"Plik nie istnieje: {file_4}")
        return

    try:
        data_states = load_data(file_4) # plik zawiera 4 stany wlasne to n max 4
        plot_n_states(data_states, out_2, 2)
        plot_n_states(data_states, out_4, 4)
    except Exception as e:
        print(f"Błąd: {e}")
        raise
 
 
if __name__ == "__main__":
    main()