#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
 
# ===================== STYL WYKRESÓW =====================
rc("text.latex", preamble=r"\usepackage{lmodern} \usepackage{physics}")
 
plt.rcParams.update({
    "font.size": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "axes.labelsize": 22,
    "axes.linewidth": 1.2,
    "lines.linewidth": 2.5,
    "figure.figsize": (10, 8),
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "text.usetex": True,
    "font.family": "Computer Modern Roman",
})
 
# ===================== STAŁE FIZYCZNE =====================
E_h    = 27.211   # [eV]
wy_meV = 200.0    # [meV] 
 
# ===================== ENERGIE ANALITYCZNE =====================
def analytical_energies(hbwx_meV, n_states=10):
    wx_eV = hbwx_meV * 1e-3
    wy_eV = wy_meV   * 1e-3
 
    levels = []
    for nx in range(20):
        for ny in range(20):
            E = wx_eV * (nx + 0.5) + wy_eV * (ny + 0.5)
            levels.append(E)
 
    levels.sort()
    return levels[:n_states]
 
 
# ===================== WCZYTYWANIE DANYCH =====================
def load_energy_data(filename):
    data = np.loadtxt(filename, comments='#')
    hbwx  = data[:, 0]
    E_num = data[:, 1:]
    return hbwx, E_num
 
 
# ===================== RYSOWANIE WYKRESU =====================
def plot_energies(hbwx, E_num, output_file):
    n_states = E_num.shape[1]
    colors   = plt.cm.tab10(np.linspace(0, 1, n_states))
 
    plt.figure()
 
    for s in range(n_states):
        plt.plot(hbwx, E_num[:, s], color=colors[s],
                 label=rf"$E_{{{s}}}$ num.")
 
    for s in range(min(n_states, 6)):
        E_anal = [analytical_energies(wx, n_states)[s] for wx in hbwx]
        plt.plot(hbwx, E_anal, color=colors[s], linestyle='--', alpha=0.7,
                 label=rf"$E_{{{s}}}$ anal." if s == 0 else None)
 
    plt.xlabel(r"$\hbar\omega_x\;[\mathrm{meV}]$")
    plt.ylabel(r"$E\;[\mathrm{eV}]$")
    plt.title(r"Energie własne w funkcji $\hbar\omega_x$, $\hbar\omega_y = 200\;\mathrm{meV}$")
    plt.legend(fontsize=11, ncol=2, loc='upper left')
    plt.grid(True, linestyle="--", alpha=0.3)
 
    plt.savefig(output_file)
    plt.close()
 
 
# ===================== MAIN =====================
def main():
    if len(sys.argv) != 3:
        print("Użycie: python3 plot_energies.py <plik_danych> <plik_wyjściowy>")
        return
 
    data_file   = Path(sys.argv[1])
    output_file = Path(sys.argv[2])
 
    if not data_file.exists():
        print("Plik nie istnieje:", data_file)
        return
 
    try:
        hbwx, E_num = load_energy_data(data_file)
        plot_energies(hbwx, E_num, output_file)
        print("Zapisano wykres:", output_file)
    except Exception as e:
        print("Błąd:", e)
 
 
if __name__ == "__main__":
    main()