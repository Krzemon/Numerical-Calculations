#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Import funkcji z zadania 1
sys.path.insert(0, str(Path(__file__).parent))
from plot_zad1 import load_data, plot_dispersion_base

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


# ===================== Mody un(y) - wszystkie na jednym wykresie =====================
def plot_modes(prefix, E_label, output_file):
    colors = ['steelblue', 'tomato', 'seagreen', 'darkorange']
    
    n_modes = 0
    for im in range(10):
        if Path(prefix + str(im) + ".dat").exists():
            n_modes += 1
        else:
            break
    
    if n_modes == 0:
        print(f"Brak plikow modow dla {prefix}")
        return
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    for im in range(n_modes):
        fname = Path(prefix + str(im) + ".dat")
        data = load_data(fname)
        y, re, im_ = data[:, 0], data[:, 1], data[:, 2]
        
        ax.plot(y, re, color=colors[im % len(colors)], 
                label=rf"$\mathrm{{Re}}\,u_{{-,{im}}}$", linewidth=2.5)
        ax.plot(y, im_, color=colors[im % len(colors)], 
                linestyle='--', label=rf"$\mathrm{{Im}}\,u_{{-,{im}}}$", linewidth=2.5)
    
    ax.set_xlabel(r"$y\;[\mathrm{nm}]$")
    ax.set_ylabel(r"$u_{-,n}(y)$")
    # ax.set_title(rf"Czesc rzeczywista i urojona funkcji wlasnych dla $E={E_label}$\,eV")
    ax.legend(fontsize=12, loc='best', ncol=2)
    ax.grid(True, linestyle="--", alpha=0.3)
    
    fig.tight_layout()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_file)
    plt.close(fig)
    print(f"Zapisano: {output_file}")


# ===================== Dyspersja z punktami dla jednej energii =====================
def plot_dispersion_with_modes(disp_data, modes_file, E_label, output_file, ylim=None):
    """Rysuje relację dyspersji z naniesionymi punktami modów"""
    kx = disp_data[:, 0]
    E = disp_data[:, 1:]

    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Używamy funkcji z zadania 1 do narysowania dyspersji
    plot_dispersion_base(ax, kx, E)

    E_val = float(E_label)
    ax.axhline(E_val, color='gray', linestyle='--', linewidth=1.0, alpha=0.6)
    
    if Path(modes_file).exists():
        modes_data = load_data(modes_file)
        if modes_data.ndim == 2 and len(modes_data) > 0:
            kx_modes = modes_data[:, 0]
            E_modes = modes_data[:, 1]
            
            # Rozdzielamy mody na dodatnie (k>0) i ujemne (k<0)
            mask_pos = kx_modes > 0
            mask_neg = kx_modes < 0
            
            # Mody dodatnie - czerwone
            if np.any(mask_pos):
                ax.plot(kx_modes[mask_pos], E_modes[mask_pos], 'ro', markersize=8, 
                        label=rf'Mody $k>0$', zorder=5)
            
            # Mody ujemne - niebieskie
            if np.any(mask_neg):
                ax.plot(kx_modes[mask_neg], E_modes[mask_neg], 'bo', markersize=8, 
                        label=rf'Mody $k<0$', zorder=5)
            
            ax.legend(fontsize=14, loc='best')

    ax.set_xlabel(r"$k_x\;[\mathrm{nm}^{-1}]$")
    ax.set_ylabel(r"$E\;[\mathrm{eV}]$")
    
    if ylim:
        ax.set_ylim(ylim)
        # ax.set_title(rf"Relacja dyspersji $E(k_x)$ z modami dla $E={E_label}$\,eV (zblizenie)")
    # else:
        # ax.set_title(rf"Relacja dyspersji $E(k_x)$ z modami dla $E={E_label}$\,eV")
    
    ax.grid(True, linestyle="--", alpha=0.3)

    fig.tight_layout()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_file)
    plt.close(fig)
    print(f"Zapisano: {output_file}")


# ===================== MAIN =====================
def main():
    if len(sys.argv) != 6:
        print("Uzycie: python3 plot_zad2.py <prefix_u> <modes.dat> <disp.dat> <E_label> <output_prefix>")
        sys.exit(1)

    prefix_u = sys.argv[1]
    modes_file = sys.argv[2]
    disp_file = sys.argv[3]
    E_label = sys.argv[4]
    out_prefix = sys.argv[5]

    disp_data = load_data(disp_file)

    # Rysujemy mody un(y)
    plot_modes(prefix_u, E_label, Path(out_prefix + "_modes.png"))
    
    # Rysujemy dyspersję z punktami - pełny zakres
    plot_dispersion_with_modes(disp_data, modes_file, E_label, 
                               Path(out_prefix + "_disp_full.png"))
    
    # Rysujemy dyspersję z punktami - zblíżenie (0-5 eV)
    plot_dispersion_with_modes(disp_data, modes_file, E_label, 
                               Path(out_prefix + "_disp_zoom.png"), ylim=(0, 5))


if __name__ == "__main__":
    main()
