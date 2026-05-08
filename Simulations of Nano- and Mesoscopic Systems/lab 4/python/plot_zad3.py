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


def load_data(filename):
    data = np.loadtxt(filename, comments='#')
    return data[~np.isnan(data).any(axis=1)]


# ===================== T(E) =====================
def plot_TE(data, output_file):
    E = data[:, 0]
    T = data[:, 1]
    R = data[:, 2]

    fig, ax = plt.subplots(figsize=(10, 7))

    ax.plot(E, T, '-', color='steelblue', label=r"$T$", linewidth=2.5)
    ax.plot(E, R, '-', color='tomato', label=r"$R$", linewidth=2.5)
    ax.plot(E, T + R, '--', color='black', label=r"$T + R$", linewidth=2, alpha=0.7)

    ax.set_xlabel(r"$E\;[\mathrm{eV}]$")
    ax.set_ylabel(r"Współczynnik")
    ax.set_title(r"Dioda RTD -- transmisja przez podwójną barierę")
    ax.legend(fontsize=14, loc='best')
    ax.set_xlim(E[0], E[-1])
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, linestyle="--", alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_file)
    plt.close(fig)
    print(f"Zapisano: {output_file}")


# ===================== I-V =====================
def plot_IV(data, output_file):
    V = data[:, 0]
    j_He = data[:, 1]
    j_N = data[:, 2]

    fig, ax = plt.subplots(figsize=(10, 7))

    ax.plot(V, j_He, '-', color='steelblue', linewidth=2.5, label=r"$T = 4.2\,\mathrm{K}$")
    ax.plot(V, j_N, '-', color='tomato', linewidth=2.5, label=r"$T = 77\,\mathrm{K}$")

    ax.set_xlabel(r"$V_{\mathrm{bias}}\;[\mathrm{mV}]$")
    ax.set_ylabel(r"$j\;[\mathrm{j.at.}]$")
    ax.set_title(r"Charakterystyka prądowo-napięciowa diody RTD")
    ax.legend(fontsize=14, loc='best')
    ax.grid(True, linestyle="--", alpha=0.3)
    ax.set_xlim(V[0], V[-1])

    fig.tight_layout()
    fig.savefig(output_file)
    plt.close(fig)
    print(f"Zapisano: {output_file}")


# ===================== MAIN =====================
def main():
    if len(sys.argv) != 5:
        print("Użycie: python3 plot_zad3.py "
              "<TE.dat> <IV.dat> <out_TE.png> <out_IV.png>")
        sys.exit(1)

    te_file  = Path(sys.argv[1])
    iv_file  = Path(sys.argv[2])
    out_te   = Path(sys.argv[3])
    out_iv   = Path(sys.argv[4])

    for f in [te_file, iv_file]:
        if not f.exists():
            print(f"Plik nie istnieje: {f}")
            sys.exit(1)

    plot_TE(load_data(te_file), out_te)
    plot_IV(load_data(iv_file), out_iv)


if __name__ == "__main__":
    main()
