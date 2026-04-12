#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy.interpolate import griddata
 
# ===================== STYL WYKRESOW =====================
rc("text.latex", preamble=r"\usepackage{lmodern} \usepackage{physics}")
 
plt.rcParams.update({
    "font.size": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "axes.labelsize": 22,
    "axes.linewidth": 1.2,
    "lines.linewidth": 2.5,
    "figure.figsize": (7, 6),
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "text.usetex": True,
    "font.family": "Computer Modern Roman",
})
 
 
# ===================== WCZYTYWANIE DANYCH =====================
def load_grid_data(filename):
    """
    Wczytuje dane z pliku tekstowego.
 
    Pierwszy wiersz komentarza zawiera "units=nm" lub "units=au",
    co decyduje o etykietach osi.
 
    Formaty plikow:
    - basis_k*.dat:  3 kolumny  x  y  phi
    - WF_state*.dat: 4 kolumny  x  y  psi  psi^2  (rysujemy kolumne 4)
 
    Zwraca X, Y, Z jako siatki 2D oraz string z jednostka ("nm" lub "au").
    """
    name  = Path(filename).stem
    units = "au"  # domyslnie
 
    data = []
    with open(filename) as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith('#'):
                if "units=nm" in stripped:
                    units = "nm"
                elif "units=au" in stripped:
                    units = "au"
                continue
            if not stripped:
                continue
            data.append([float(v) for v in stripped.split()])
 
    data = np.array(data)
    x = data[:, 0]
    y = data[:, 1]
 
    # WF_state* ma 4 kolumny: x y psi psi^2 -- rysujemy psi^2
    if name.startswith("WF_state") and data.shape[1] >= 4:
        z = data[:, 3]
    else:
        z = data[:, 2]
 
    ux = np.unique(np.round(x, 10))
    uy = np.unique(np.round(y, 10))
    nx, ny = len(ux), len(uy)
 
    if nx * ny == len(x):
        X = x.reshape((nx, ny))
        Y = y.reshape((nx, ny))
        Z = z.reshape((nx, ny))
    else:
        xi = np.linspace(x.min(), x.max(), 100)
        yi = np.linspace(y.min(), y.max(), 100)
        X, Y = np.meshgrid(xi, yi)
        Z = griddata((x, y), z, (X, Y), method='cubic')
 
    return X, Y, Z, units
 
 
# ===================== TYTUL I ETYKIETY =====================
def make_title(data_file: Path) -> str:
    name = data_file.stem
    if name.startswith("basis_k"):
        k = name.split("basis_k")[1]
        return rf"$\phi_{{k={k}}}(x,y)$"
    elif name.startswith("WF_state") and name.endswith("_wy"):
        s = name.replace("WF_state", "").replace("_wy", "")
        return rf"$|\Psi_{{{s}}}|^2$, $\hbar\omega_y=400\,\mathrm{{meV}}$"
    elif name.startswith("WF_state"):
        s = name.replace("WF_state", "")
        return rf"$|\Psi_{{{s}}}|^2$"
    else:
        return name
 
 
def axis_label(units: str) -> tuple:
    """Zwraca (xlabel, ylabel) dla danej jednostki."""
    if units == "nm":
        return r"$x\;[\mathrm{nm}]$", r"$y\;[\mathrm{nm}]$"
    else:
        return r"$x\;[a_0]$", r"$y\;[a_0]$"
 
 
# ===================== RYSOWANIE WYKRESU =====================
def plot_contour(X, Y, Z, units, output_file, title):
    plt.figure()
 
    contour = plt.contourf(X, Y, Z, levels=60, cmap="viridis") # cmap="viridis") , cmap="RdBu_r"
 
    cbar = plt.colorbar(contour)
    cbar.set_label(r"$\phi(x,y)$", fontsize=18)
 
    xlabel, ylabel = axis_label(units)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
 
    plt.gca().set_aspect("equal", adjustable="box")
    plt.grid(True, linestyle="--", alpha=0.3)
 
    plt.savefig(output_file)
    plt.close()
 
 
# ===================== MAIN =====================
def main():
    if len(sys.argv) < 3:
        print("Uzycie: python3 plot_contour.py <plik_danych> <plik_wyjsciowy>")
        return
 
    data_file   = Path(sys.argv[1])
    output_file = Path(sys.argv[2])
    title       = make_title(data_file)
 
    if not data_file.exists():
        print("Plik nie istnieje:", data_file)
        return
 
    try:
        X, Y, Z, units = load_grid_data(data_file)
        plot_contour(X, Y, Z, units, output_file, title)
        print("Zapisano wykres:", output_file)
    except Exception as e:
        print("Blad:", e)
 
 
if __name__ == "__main__":
    main()