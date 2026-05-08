#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

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
    data = np.loadtxt(filename, comments='#')
    return data[~np.isnan(data).any(axis=1)]

# ===================== RYSOWANIE =====================
def plot_zad1(data, output_dir):
    E = data[:, 0]
    T_num = data[:, 1]
    R_num = data[:, 2]
    T_ana = data[:, 3]
    R_ana = data[:, 4]

    # Wykres 1: T, R numeryczne
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(E, T_num, '-', color='steelblue', label=r"$T$", linewidth=2.5)
    ax.plot(E, R_num, '-', color='tomato', label=r"$R$", linewidth=2.5)
    ax.plot(E, T_num + R_num, '--', color='black', label=r"$T + R$", linewidth=2, alpha=0.7)
    ax.set_xlabel(r"$E\;[\mathrm{eV}]$")
    ax.set_ylabel(r"Współczynnik")
    ax.set_title(r"Stała masa $m^* = m^*_{\mathrm{GaAs}}$ -- numeryczna metoda macierzy transferu")
    ax.legend(fontsize=14, loc='best')
    ax.set_xlim(E[0], E[-1])
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, linestyle="--", alpha=0.3)
    fig.tight_layout()
    out1 = output_dir / "TR_numerical.png"
    fig.savefig(out1)
    plt.close(fig)
    print(f"Zapisano: {out1}")

    # Wykres 2: T, R analityczne
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(E, T_ana, '-', color='steelblue', label=r"$T$", linewidth=2.5)
    ax.plot(E, R_ana, '-', color='tomato', label=r"$R$", linewidth=2.5)
    ax.plot(E, T_ana + R_ana, '--', color='black', label=r"$T + R$", linewidth=2, alpha=0.7)
    ax.set_xlabel(r"$E\;[\mathrm{eV}]$")
    ax.set_ylabel(r"Współczynnik")
    ax.set_title(r"Stała masa $m^* = m^*_{\mathrm{GaAs}}$ -- rozwiązanie analityczne")
    ax.legend(fontsize=14, loc='best')
    ax.set_xlim(E[0], E[-1])
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, linestyle="--", alpha=0.3)
    fig.tight_layout()
    out2 = output_dir / "TR_analytical.png"
    fig.savefig(out2)
    plt.close(fig)
    print(f"Zapisano: {out2}")

    # Wykres 3: Różnica
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(E, np.abs(T_num - T_ana), '-', color='steelblue', label=r"$|T_{\mathrm{num}} - T_{\mathrm{ana}}|$", linewidth=2.5)
    ax.plot(E, np.abs(R_num - R_ana), '-', color='tomato', label=r"$|R_{\mathrm{num}} - R_{\mathrm{ana}}|$", linewidth=2.5)
    ax.set_xlabel(r"$E\;[\mathrm{eV}]$")
    ax.set_ylabel(r"Błąd bezwzględny")
    ax.set_title(r"Różnica między metodą numeryczną a analityczną")
    ax.legend(fontsize=14, loc='best')
    ax.set_xlim(E[0], E[-1])
    ax.set_yscale('log')
    ax.grid(True, linestyle="--", alpha=0.3)
    fig.tight_layout()
    out3 = output_dir / "TR_difference.png"
    fig.savefig(out3)
    plt.close(fig)
    print(f"Zapisano: {out3}")

def plot_zad2(data, output_dir):
    E = data[:, 0]
    T_var = data[:, 1]
    R_var = data[:, 2]

    # Wykres 1: T, R dla zmiennej masy numerycznie
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(E, T_var, '-', color='steelblue', label=r"$T$", linewidth=2.5)
    ax.plot(E, R_var, '-', color='tomato', label=r"$R$", linewidth=2.5)
    ax.plot(E, T_var + R_var, '--', color='black', label=r"$T + R$", linewidth=2, alpha=0.7)
    ax.set_xlabel(r"$E\;[\mathrm{eV}]$")
    ax.set_ylabel(r"Współczynnik")
    ax.set_title(r"Zmienna masa $m^*(z)$ -- metoda macierzy transferu")
    ax.legend(fontsize=14, loc='best')
    ax.set_xlim(E[0], E[-1])
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, linestyle="--", alpha=0.3)
    fig.tight_layout()
    out1 = output_dir / "TR_var_mass.png"
    fig.savefig(out1)
    plt.close(fig)
    print(f"Zapisano: {out1}")

# ===================== MAIN =====================
def main():
    if len(sys.argv) != 5:
        print("Użycie: python3 plot_zad1_2.py <data_zad1.dat> <data_zad2.dat> <out_dir_zad1> <out_dir_zad2>")
        sys.exit(1)

    data_file_zad1 = Path(sys.argv[1])
    data_file_zad2 = Path(sys.argv[2])
    out_dir_zad1 = Path(sys.argv[3])
    out_dir_zad2 = Path(sys.argv[4])

    if not data_file_zad1.exists():
        print(f"Plik nie istnieje: {data_file_zad1}")
        sys.exit(1)
    
    if not data_file_zad2.exists():
        print(f"Plik nie istnieje: {data_file_zad2}")
        sys.exit(1)

    out_dir_zad1.mkdir(parents=True, exist_ok=True)
    out_dir_zad2.mkdir(parents=True, exist_ok=True)

    data_zad1 = load_data(data_file_zad1)
    data_zad2 = load_data(data_file_zad2)
    
    plot_zad1(data_zad1, out_dir_zad1)
    plot_zad2(data_zad2, out_dir_zad2)

if __name__ == "__main__":
    main()
