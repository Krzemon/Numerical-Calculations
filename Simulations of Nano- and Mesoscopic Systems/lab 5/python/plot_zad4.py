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


# ===================== T(Vgates) - 3 osobne wykresy =====================
def plot_conductance(data, output_prefix):
    Vg = data[:, 0]
    T = data[:, 1]
    R = data[:, 2]
    TR = data[:, 3]

    # (a) T(Vgates)
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(Vg, T, color='steelblue', linewidth=2.5)
    
    # Linie poziome dla kwantowania (n = 0, 1, 2, ...)
    max_T = int(np.ceil(np.max(T)))
    for n in range(max_T + 1):
        ax.axhline(n, color='gray', linestyle=':', linewidth=1.2, alpha=0.6)
        ax.text(Vg[-1] + 0.01, n, f'$n={n}$', fontsize=12, va='center')
    
    ax.set_xlabel(r"$V_\mathrm{gates}\;[\mathrm{eV}]$")
    ax.set_ylabel(r"Transmisja $T$ [$2e^2/h$]")
    # ax.set_title(r"(a) Transmisja $T(V_\mathrm{gates})$ -- kwantowanie przewodności")
    ax.grid(True, linestyle="--", alpha=0.3)
    ax.set_xlim(Vg[0], Vg[-1])
    ax.set_ylim(-0.2, max_T + 0.5)
    fig.tight_layout()
    out_T = Path(str(output_prefix) + "_T.png")
    out_T.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_T)
    plt.close(fig)
    print(f"Zapisano: {out_T}")

    # (b) R(Vgates)
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(Vg, R, color='tomato', linewidth=2.5)
    ax.set_xlabel(r"$V_\mathrm{gates}\;[\mathrm{eV}]$")
    ax.set_ylabel(r"Odbicie $R$ [$2e^2/h$]")
    # ax.set_title(r"(b) Odbicie $R(V_\mathrm{gates})$")
    ax.grid(True, linestyle="--", alpha=0.3)
    ax.set_xlim(Vg[0], Vg[-1])
    fig.tight_layout()
    out_R = Path(str(output_prefix) + "_R.png")
    fig.savefig(out_R)
    plt.close(fig)
    print(f"Zapisano: {out_R}")

    # (c) T+R(Vgates)
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(Vg, TR, color='black', linewidth=2.5)
    
    # Linie poziome dla całkowitej liczby modów
    max_TR = int(np.ceil(np.max(TR)))
    for n in range(max_TR + 1):
        ax.axhline(n, color='gray', linestyle=':', linewidth=1.2, alpha=0.6)
    
    ax.set_xlabel(r"$V_\mathrm{gates}\;[\mathrm{eV}]$")
    ax.set_ylabel(r"$T+R$ [liczba modów]")
    # ax.set_title(r"(c) Suma $T(V_\mathrm{gates}) + R(V_\mathrm{gates})$ -- zachowanie liczby modów")
    ax.grid(True, linestyle="--", alpha=0.3)
    ax.set_xlim(Vg[0], Vg[-1])
    fig.tight_layout()
    out_TR = Path(str(output_prefix) + "_TR.png")
    fig.savefig(out_TR)
    plt.close(fig)
    print(f"Zapisano: {out_TR}")


# ===================== Mapa V(x,y) =====================
def plot_potential(data, output_file):
    x, y, V = data[:, 0], data[:, 1], data[:, 2]

    xs = np.unique(x)
    ys = np.unique(y)
    Z = V.reshape(len(xs), len(ys)).T

    # Interpolacja dla gladszego wykresu - zakomentowane
    # Wartosci x,y pochodza z dyskretyzacji siatki w QPC.cpp (make_QTBM)
    # from scipy.interpolate import RectBivariateSpline
    # interp = RectBivariateSpline(ys, xs, Z)
    # xs_fine = np.linspace(xs.min(), xs.max(), 300)
    # ys_fine = np.linspace(ys.min(), ys.max(), 300)
    # Z_fine = interp(ys_fine, xs_fine)

    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.pcolormesh(xs, ys, Z, shading='auto', cmap='inferno')
    cbar = fig.colorbar(im, ax=ax, label=r"$V\;[\mathrm{meV}]$")
    ax.set_xlabel(r"$x\;[\mathrm{nm}]$")
    ax.set_ylabel(r"$y\;[\mathrm{nm}]$")
    # ax.set_title(r"Potencjal $V_\mathrm{QPC}(x,y)$ dla $V_\mathrm{gates}=-1\,\mathrm{eV}$")
    ax.set_aspect('equal')

    fig.tight_layout()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_file)
    plt.close(fig)
    print(f"Zapisano: {output_file}")


# ===================== MAIN =====================
def main():
    if len(sys.argv) != 5:
        print("Użycie: python3 plot_zad4.py <conductance.dat> <potential.dat> "
              "<out_conductance_prefix> <out_potential.png>")
        sys.exit(1)

    cond_file = Path(sys.argv[1])
    pot_file = Path(sys.argv[2])
    out_cond_prefix = sys.argv[3]
    out_pot = Path(sys.argv[4])

    for f in [cond_file, pot_file]:
        if not f.exists():
            print(f"Plik nie istnieje: {f}")
            sys.exit(1)

    plot_conductance(load_data(cond_file), out_cond_prefix)
    plot_potential(load_data(pot_file), out_pot)


if __name__ == "__main__":
    main()
