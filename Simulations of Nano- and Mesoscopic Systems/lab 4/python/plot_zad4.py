#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap

# ===================== STYL =====================
plt.rcParams["text.usetex"] = True
plt.rcParams["text.latex.preamble"] = r"""
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{lmodern}
\usepackage{amsmath}
"""
plt.rcParams.update({
    "font.size": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "axes.labelsize": 18,
    "axes.linewidth": 1.2,
    "lines.linewidth": 2.5,
    "figure.figsize": (10, 6),
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "font.family": "serif",
})


def load_data(filename):
    return np.loadtxt(filename, comments='#')


COLORS = ['steelblue', 'tomato', 'seagreen', 'darkorange', 'purple']
STYLES = ['-', '--', '-.', ':', (0, (3, 1, 1, 1))]


# ===================== E_n(x) =====================
def plot_En(data, output_file):
    x_nm    = data[:, 0]
    n_chan   = data.shape[1] - 1

    fig, ax = plt.subplots(figsize=(11, 7))

    for n in range(n_chan):
        En_meV = data[:, n + 1]
        color  = COLORS[n % len(COLORS)]
        ls     = STYLES[n % len(STYLES)]
        ax.plot(x_nm, En_meV, linestyle=ls, color=color,
                label=r"$E_{" + str(n + 1) + r"}(x)$")

    ax.set_xlabel(r"$x\;[\mathrm{nm}]$")
    ax.set_ylabel(r"$E_n(x)\;[\mathrm{meV}]$")
    ax.set_title(
        r"Efektywne potencjały QPC w przybliżeniu adiabatycznym"
        "\n"
        r"($V_g = 4\,\mathrm{eV}$, $W = 50\,\mathrm{nm}$, $L = 100\,\mathrm{nm}$)"
    )
    ax.legend(fontsize=13, ncol=2)
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.set_xlim(x_nm[0], x_nm[-1])

    fig.tight_layout()
    fig.savefig(output_file)
    plt.close(fig)
    print(f"Zapisano: {output_file}")


# ===================== G(E) =====================
def plot_GE(data, output_file):
    E_meV = data[:, 0]
    G     = data[:, 1]

    n_max = int(G.max()) + 2

    fig, ax = plt.subplots(figsize=(10, 7))

    ax.plot(E_meV, G, color='steelblue', label=r"$G(E)$")

    for n in range(1, n_max):
        ax.axhline(n, color='gray', linestyle=':', linewidth=1.2, alpha=0.6)
        ax.text(E_meV[-1] * 1.01, n, rf"${n} \times 2e^2/h$",
                va='center', fontsize=11, color='gray')

    ax.set_xlabel(r"$E\;[\mathrm{meV}]$")
    ax.set_ylabel(r"$G\;[2e^2/h]$")
    ax.set_title(
        r"Kwantyzacja konduktancji QPC"
        "\n"
        r"(formuła Landauera, przybliżenie adiabatyczne)"
    )
    ax.legend(fontsize=13)
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.set_xlim(E_meV[0], E_meV[-1] * 1.01)
    ax.set_ylim(-0.1, n_max - 0.5)

    fig.tight_layout()
    fig.savefig(output_file)
    plt.close(fig)
    print(f"Zapisano: {output_file}")


# ===================== G(Vg) =====================
def plot_GVg(data, output_file):
    Vg_meV = data[:, 0]
    G1    = data[:, 1]
    G2    = data[:, 2]

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.plot(Vg_meV / 1000.0, G1, 'b-',  label=r"$E_F = 50\,\mathrm{meV}$")
    ax.plot(Vg_meV / 1000.0, G2, 'r--', label=r"$E_F = 100\,\mathrm{meV}$")

    # Poziomy kwantyzacji
    n_max = int(max(G1.max(), G2.max())) + 2
    for n in range(1, n_max):
        ax.axhline(n, color='gray', linestyle=':', linewidth=1.0, alpha=0.5)

    ax.set_xlabel(r"$V_g\;[\mathrm{V}]$")
    ax.set_ylabel(r"$G\;[2e^2/h]$")
    ax.set_title(
        r"Konduktancja QPC w funkcji napięcia bramki"
        "\n"
        r"(formuła Landauera)"
    )
    ax.legend(fontsize=13)
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.set_xlim(Vg_meV[0] / 1000.0, Vg_meV[-1] / 1000.0 * 1.01)

    fig.tight_layout()
    fig.savefig(output_file)
    plt.close(fig)
    print(f"Zapisano: {output_file}")


# ===================== MAIN =====================
def main():
    if len(sys.argv) != 7:
        print("Użycie: python3 plot_zad4.py "
              "<En.dat> <GE.dat> <GVg.dat> "
              "<out_En.png> <out_GE.png> <out_GVg.png>")
        sys.exit(1)

    en_file  = Path(sys.argv[1])
    ge_file  = Path(sys.argv[2])
    gvg_file = Path(sys.argv[3])
    out_en   = Path(sys.argv[4])
    out_ge   = Path(sys.argv[5])
    out_gvg  = Path(sys.argv[6])

    for f in [en_file, ge_file, gvg_file]:
        if not f.exists():
            print(f"Plik nie istnieje: {f}")
            sys.exit(1)

    plot_En (load_data(en_file),  out_en)
    plot_GE (load_data(ge_file),  out_ge)
    plot_GVg(load_data(gvg_file), out_gvg)


if __name__ == "__main__":
    main()
