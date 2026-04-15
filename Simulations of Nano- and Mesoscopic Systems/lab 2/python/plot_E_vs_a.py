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
    "figure.figsize": (10, 7),
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "text.usetex": True,
    "font.family": "Computer Modern Roman",
})
 
 
# ===================== WCZYTYWANIE DANYCH =====================
def load_data(filename):
    """Wczytuje dwie kolumny: a [nm], E [eV]."""
    data = np.loadtxt(filename, comments='#')
    return data[:, 0], data[:, 1]
 
 
# ===================== RYSOWANIE WYKRESU =====================
def plot_energies(a_IT, E_IT, a_HF, E_HF, output_file):
    fig, ax = plt.subplots()
 
    ax.plot(a_IT, E_IT * 1e3, 'o-',  color='steelblue', label=r"Czas urojony (dokladne)")
    ax.plot(a_HF, E_HF * 1e3, 's--', color='tomato',    label=r"Hartree-Fock (pole srednie)")
 
    ax.set_xlabel(r"$a\;[\mathrm{nm}]$")
    ax.set_ylabel(r"$E\;[\mathrm{meV}]$")
    ax.set_title(r"Energia stanu podstawowego 2 elektronow od rozmiaru kropki")
    ax.legend(fontsize=14)
    ax.grid(True, linestyle="--", alpha=0.3)
 
    fig.savefig(output_file)
    plt.close(fig)
    print(f"Zapisano wykres: {output_file}")
 
 
# ===================== MAIN =====================
def main():
    if len(sys.argv) != 4:
        print("Uzycie: python3 plot_E_vs_a.py <plik_IT> <plik_HF> <plik_wyjsciowy>")
        return
 
    file_IT     = Path(sys.argv[1])
    file_HF     = Path(sys.argv[2])
    output_file = Path(sys.argv[3])
 
    for f in (file_IT, file_HF):
        if not f.exists():
            print(f"Plik nie istnieje: {f}")
            return
 
    try:
        a_IT, E_IT = load_data(file_IT)
        a_HF, E_HF = load_data(file_HF)
        plot_energies(a_IT, E_IT, a_HF, E_HF, output_file)
    except Exception as e:
        print(f"Blad: {e}")
        raise
 
 
if __name__ == "__main__":
    main()
    