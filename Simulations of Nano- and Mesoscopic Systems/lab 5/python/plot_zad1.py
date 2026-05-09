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
    if data.ndim == 1:
        return data
    return data[~np.isnan(data).any(axis=1)]


def plot_dispersion_base(ax, kx, E):
    """Rysuje relację dyspersji na podanym axes"""
    for i in range(E.shape[1]):
        ax.plot(kx, E[:, i], color='steelblue', linewidth=1.5)


# ===================== E(kx) bez punktów =====================
def plot_dispersion(data, output_file):
    kx = data[:, 0]
    E = data[:, 1:]

    fig, ax = plt.subplots(figsize=(10, 7))
    plot_dispersion_base(ax, kx, E)

    ax.set_xlabel(r"$k_x\;[\mathrm{nm}^{-1}]$")
    ax.set_ylabel(r"$E\;[\mathrm{eV}]$")
    # ax.set_title(r"Relacja dyspersji $E(k_x)$ -- jednorodny kanał")
    ax.grid(True, linestyle="--", alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_file)
    plt.close(fig)
    print(f"Zapisano: {output_file}")


# ===================== E(kx) z punktami dla jednej energii =====================
def plot_dispersion_with_modes(disp_data, modes_file, E_label, output_file):
    kx = disp_data[:, 0]
    E = disp_data[:, 1:]

    fig, ax = plt.subplots(figsize=(10, 7))
    plot_dispersion_base(ax, kx, E)

    E_val = float(E_label)
    ax.axhline(E_val, color='gray', linestyle='--', linewidth=1.0, alpha=0.6)
    
    if Path(modes_file).exists():
        modes_data = load_data(modes_file)
        if modes_data.ndim == 2 and len(modes_data) > 0:
            kx_modes = modes_data[:, 0]
            
            # Mody propagujace w prawo (u-, dodatnie kx)
            E_pos = np.full(len(kx_modes), E_val)
            ax.plot(kx_modes, E_pos, 'ro', markersize=6, 
                    label=rf'$u_-$ (prawo, $+k_x$)', zorder=5)
            
            # Mody propagujace w lewo (u+, ujemne kx)
            kx_neg = -kx_modes[::-1]
            E_neg = np.full(len(kx_neg), E_val)
            ax.plot(kx_neg, E_neg, 'bo', markersize=6, 
                    label=rf'$u_+$ (lewo, $-k_x$)', zorder=5)
            
            ax.legend(fontsize=14, loc='best')

    ax.set_xlabel(r"$k_x\;[\mathrm{nm}^{-1}]$")
    ax.set_ylabel(r"$E\;[\mathrm{eV}]$")
    ax.set_title(rf"Relacja dyspersji $E(k_x)$ z modami dla $E={E_label}$\,eV")
    ax.grid(True, linestyle="--", alpha=0.3)

    fig.tight_layout()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_file)
    plt.close(fig)
    print(f"Zapisano: {output_file}")


# ===================== MAIN =====================
def main():
    if len(sys.argv) < 3:
        print("Użycie:")
        print("  python3 plot_zad1.py <dispersion.dat> <output.png>")
        print("  python3 plot_zad1.py <dispersion.dat> <modes.dat> <E_label> <output.png>")
        sys.exit(1)

    data_file = Path(sys.argv[1])
    
    if not data_file.exists():
        print(f"Plik nie istnieje: {data_file}")
        sys.exit(1)

    disp_data = load_data(data_file)

    if len(sys.argv) == 3:
        # Tylko dyspersja bez punktów
        out_file = Path(sys.argv[2])
        out_file.parent.mkdir(parents=True, exist_ok=True)
        plot_dispersion(disp_data, out_file)
    elif len(sys.argv) == 5:
        # Dyspersja z punktami
        modes_file = sys.argv[2]
        E_label = sys.argv[3]
        out_file = Path(sys.argv[4])
        plot_dispersion_with_modes(disp_data, modes_file, E_label, out_file)
    else:
        print("Błędna liczba argumentów")
        sys.exit(1)


if __name__ == "__main__":
    main()
