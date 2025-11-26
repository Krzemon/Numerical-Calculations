#!/usr/bin/env python3
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

def load_data(filename):
    data = np.loadtxt(filename, comments='#')
    return {
        "element": data[:, 0].astype(int),
        "global_node": data[:, 1].astype(int),
        "local_node": data[:, 2].astype(int),
        "x": data[:, 3],
        "y": data[:, 4],
    }

def main():
    if len(sys.argv) != 4:
        print("Uzycie: python3 plot_nodes.py <plik_danych> <plik_wykresu> <tytul>")
        return

    data_file = Path(sys.argv[1])
    out_png = Path(sys.argv[2])
    title = sys.argv[3]

    if not data_file.exists():
        print("Plik nie istnieje:", data_file)
        return

    nodes = load_data(data_file)

    element = nodes["element"]
    global_node = nodes["global_node"]
    local_node = nodes["local_node"]
    x = nodes["x"]
    y = nodes["y"]

    plt.figure(figsize=(8, 6))

    # Węzły globalne
    plt.scatter(x, y, color="blue", label="Globalne", zorder=3)

    # Rysowanie prostokątów oddzielających elementy
    max_elem = int(np.max(element))
    for idx in range(1, max_elem + 1):
        ex = x[element == idx]
        ey = y[element == idx]
        min_x, max_x = ex.min(), ex.max()
        min_y, max_y = ey.min(), ey.max()
        rect = Rectangle(
            (min_x, min_y),
            max_x - min_x,
            max_y - min_y,
            linewidth=1,
            edgecolor="gray",
            facecolor="none",
            linestyle="-",
            zorder=2
        )
        plt.gca().add_patch(rect)

    # Numery globalnych węzłów
    for xn, yn, idx in zip(x, y, global_node):
        plt.annotate(
            int(idx),
            (xn, yn),
            textcoords="offset points",
            xytext=(0, 5),
            ha="center",
            fontsize=10,
            color="blue"
        )

    # Numery lokalnych węzłów przesunięte w kierunku środka elementu
    for idx_elem in range(1, max_elem + 1):
        ex = x[element == idx_elem]
        ey = y[element == idx_elem]
        lx = x[element == idx_elem]
        ly = y[element == idx_elem]
        # Środek elementu
        center_x = np.mean(ex)
        center_y = np.mean(ey)
        # przesunięcie od węzła w stronę środka elementu o 20% 
        for xn, yn, lidx in zip(lx, ly, local_node[element == idx_elem]):
            dx = (center_x - xn) * 0.2
            dy = (center_y - yn) * 0.2
            plt.annotate(
                int(lidx),
                (xn + dx, yn + dy),
                textcoords="offset points",
                xytext=(0, 0),
                ha="center",
                va="center",
                fontsize=12,
                color="darkgreen",
                zorder=5
            )

    # Numery elementów w środku
    for idx in range(1, max_elem + 1):
        center_x = np.mean(x[element == idx])
        center_y = np.mean(y[element == idx])
        plt.annotate(
            idx,
            (center_x, center_y),
            textcoords="offset points",
            xytext=(0, -2),
            ha="center",
            fontsize=12,
            color="red"
        )

    plt.title(title)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    print("Zapisano wykres:", out_png)
    plt.close()

if __name__ == "__main__":
    main()