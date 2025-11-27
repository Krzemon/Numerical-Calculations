#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def load_data_with_header(filename):
    """Wczytuje dane i wyciąga n=... z nagłówka."""
    with open(filename, "r") as f:
        header = f.readline().strip()

    # Format: "# x y u_numerical u_analytical n=3"
    if "n=" in header:
        n = header.split("n=")[-1]
    else:
        n = "?"

    data = np.loadtxt(filename, comments="#")
    x, y, u_num, u_exact = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    return x, y, u_num, u_exact, n


def prep_image_data(x, y, u):
    """Przygotowuje dane do contourf (zamiana na siatkę 2D)."""
    xs = np.unique(x)
    ys = np.unique(y)
    X, Y = np.meshgrid(xs, ys)

    U = u.reshape(len(ys), len(xs))
    extent = [xs.min(), xs.max(), ys.min(), ys.max()]
    return U, extent


def main():
    if len(sys.argv) != 3:
        print("Uzycie: python3 plot_nodes.py <plik_danych> <plik_wykresu>")
        return

    data_file = Path(sys.argv[1])
    out_png = Path(sys.argv[2])

    if not data_file.exists():
        print("Plik nie istnieje:", data_file)
        return

    # Wczytaj dane
    x, y, u_num, u_exact, n = load_data_with_header(data_file)

    # Przygotowanie danych
    U_num, extent = prep_image_data(x, y, u_num)
    U_exact, _ = prep_image_data(x, y, u_exact)

    # ==== TWORZENIE DWÓCH WYKRESÓW ====
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))

    # Numerical
    ax = axes[0]
    c1 = ax.contourf(U_num, extent=extent, cmap="RdBu", levels=20)
    fig.colorbar(c1, ax=ax)
    ax.set_title(f"Numerical solution (n={n})")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    # Analytical
    ax = axes[1]
    c2 = ax.contourf(U_exact, extent=extent, cmap="RdBu", levels=20)
    fig.colorbar(c2, ax=ax)
    ax.set_title(f"Analytical solution (n={n})")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    print("Zapisano wykres:", out_png)

    plt.close()


if __name__ == "__main__":
    main()