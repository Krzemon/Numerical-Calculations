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
    """Kolumny: t[ns], p0, p1, psum."""
    data = np.loadtxt(filename, comments='#')
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3]
 
 
# ===================== RYSOWANIE =====================
def plot_evolution(t, p0, p1, psum, title, output_file):

    fig, ax = plt.subplots()

    ax.plot(t, p0,   '-',  color='steelblue',
            label=r"$|\langle\Psi(t)|0\rangle|^2$")
    ax.plot(t, p1,   '--', color='tomato',
            label=r"$|\langle\Psi(t)|1\rangle|^2$")
    ax.plot(t, psum, ':',  color='seagreen', linewidth=1.5,
            label=r"suma")

    ax.set_xlabel(r"$t\;[\mathrm{ns}]$")
    ax.set_ylabel("Prawdopodobieństwo")
    ax.set_title(title)

    ax.set_ylim(-0.05, 1.15)
    ax.legend(fontsize=14)
    ax.grid(True, linestyle="--", alpha=0.3)

    fig.savefig(output_file)
    plt.close(fig)
 
 
# ===================== MAIN =====================
def main():
    if len(sys.argv) != 5:
        print("Użycie: python3 plot_evolution.py " 
                        "<resonance.dat> <95percent.dat> <out_res.png> <out_95.png>")
        return
 
    file_res = Path(sys.argv[1])
    file_95  = Path(sys.argv[2])
    out_res  = Path(sys.argv[3])
    out_95   = Path(sys.argv[4])
 
    for f in (file_res, file_95):
        if not f.exists():
            print(f"Plik nie istnieje: {f}")
            return
 
    try:
        t, p0, p1, psum = load_data(file_res)
        plot_evolution(
            t, p0, p1, psum,
            r"Ewolucja przy częstości rezonansowej "
            r"$\omega = \Delta E / \hbar$",
            out_res
        )

        t, p0, p1, psum = load_data(file_95)
        plot_evolution(
            t, p0, p1, psum,
            r"Ewolucja przy $\omega = 0.95\,\Delta E / \hbar$",
            out_95
        )
    except Exception as e:
        print(f"Błąd: {e}")
        raise
 
 
if __name__ == "__main__":
    main()