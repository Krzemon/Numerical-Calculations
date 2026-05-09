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


# ===================== Mapa |psi|^2 =====================
def plot_wavefunction(data, E_label, output_file):
    x, y, psi2 = data[:, 0], data[:, 1], data[:, 2]

    xs = np.unique(x)
    ys = np.unique(y)
    Z = psi2.reshape(len(xs), len(ys)).T

    # Interpolacja - zakomentowane
    # Wartosci x,y pochodza z dyskretyzacji siatki w QPC.cpp (transmission)
    # Siatka: x_min + (i+1)*dx, y_min + (j+1)*dy dla i=0..Nx-1, j=0..Ny-1
    # from scipy.interpolate import RectBivariateSpline
    # interp = RectBivariateSpline(ys, xs, Z)
    # xs_fine = np.linspace(xs.min(), xs.max(), 300)
    # ys_fine = np.linspace(ys.min(), ys.max(), 300)
    # Z_fine = interp(ys_fine, xs_fine)

    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.pcolormesh(xs, ys, Z, shading='auto', cmap='inferno')
    cbar = fig.colorbar(im, ax=ax, label=r"$|\psi|^2$")
    ax.set_xlabel(r"$x\;[\mathrm{nm}]$")
    ax.set_ylabel(r"$y\;[\mathrm{nm}]$")
    # ax.set_title(rf"Kwadrat modulu funkcji falowej dla $E={E_label}$\,eV")

    fig.tight_layout()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)
    print(f"Zapisano: {output_file}")


# ===================== MAIN =====================
def main():
    if len(sys.argv) != 4:
        print("Użycie: python3 plot_zad3.py <wavefunction.dat> <E_label> <output.png>")
        sys.exit(1)

    wf_file = Path(sys.argv[1])
    E_label = sys.argv[2]
    out_file = Path(sys.argv[3])

    if not wf_file.exists():
        print(f"Plik nie istnieje: {wf_file}")
        sys.exit(1)

    plot_wavefunction(load_data(wf_file), E_label, out_file)


if __name__ == "__main__":
    main()
