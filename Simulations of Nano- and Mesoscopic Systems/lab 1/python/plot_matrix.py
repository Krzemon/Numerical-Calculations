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
    "axes.labelsize": 18,
    "axes.linewidth": 1.2,
    "figure.figsize": (7, 6),
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "text.usetex": True,
    "font.family": "Computer Modern Roman",
})
 
 
# ===================== WCZYTYWANIE DANYCH =====================
def load_matrix(filename):
    """
    Wczytuje macierz z pliku w formacie trojkolumnowym: i  j  M(i,j).
    Pierwszy wiersz komentarza zawiera "size=NxN".
    Zwraca macierz 2D i nazwe (np. "H" lub "S").
    """
    name   = "M"
    size   = None
 
    with open(filename) as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("# matrix="):
                # "# matrix=H  size=81x81"
                parts = stripped.split()
                name  = parts[1].split("=")[1]
                for p in parts:
                    if p.startswith("size="):
                        size = int(p.split("=")[1].split("x")[0])
                break
 
    data = np.loadtxt(filename, comments='#')
    i = data[:, 0].astype(int)
    j = data[:, 1].astype(int)
    v = data[:, 2]
 
    if size is None:
        size = max(i.max(), j.max()) + 1
 
    M = np.zeros((size, size))
    M[i, j] = v
 
    return M, name
 
 
# ===================== RYSOWANIE MAPY CIEPLA =====================
def plot_heatmap(M, name, output_file):
    plt.figure()
 
    im = plt.imshow(M, cmap="RdBu_r", aspect="equal",
                    interpolation="nearest",
                    vmin=-np.abs(M).max(), vmax=np.abs(M).max())
 
    cbar = plt.colorbar(im)
    cbar.set_label(rf"${name}_{{kl}}$", fontsize=18)
 
    plt.xlabel(r"$l$")
    plt.ylabel(r"$k$")
    plt.title(rf"Macierz $\hat{{{name}}}$")
 
    plt.savefig(output_file)
    plt.close()
    print("Zapisano wykres:", output_file)
 
 
# ===================== MAIN =====================
def main():
    if len(sys.argv) < 3:
        print("Uzycie: python3 plot_matrix.py <plik_danych> <plik_wyjsciowy>")
        return
 
    data_file   = Path(sys.argv[1])
    output_file = Path(sys.argv[2])
 
    if not data_file.exists():
        print("Plik nie istnieje:", data_file)
        return
 
    try:
        M, name = load_matrix(data_file)
        plot_heatmap(M, name, output_file)
    except Exception as e:
        print("Blad:", e)
 
 
if __name__ == "__main__":
    main()