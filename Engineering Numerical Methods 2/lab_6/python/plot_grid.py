#!/usr/bin/env python3
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

def load_sections(filename):
    nodes = []
    elements = []
    nx = ny = None

    section = None
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if line == "":
                section = None
                continue

            if line.startswith("#"):
                if "index" in line:
                    section = "nodes"
                    # pobranie nx i ny
                    parts = line.split()
                    for p in parts:
                        if p.startswith("nx="):
                            nx = int(p.split("=")[1])
                        if p.startswith("ny="):
                            ny = int(p.split("=")[1])
                elif "element" in line:
                    section = "elements"
                continue

            parts = line.split()
            if section == "nodes":
                idx, x, y = parts
                nodes.append((int(idx), float(x), float(y)))
            elif section == "elements":
                elem_id, i, j, k = parts
                elements.append((int(elem_id), int(i), int(j), int(k)))

    return np.array(nodes), np.array(elements), nx, ny

def main():
    if len(sys.argv) != 4:
        print("Użycie: python3 plot_zad1_grid.py <plik_danych> <plik_png> <tytuł>")
        return

    data_file = Path(sys.argv[1])
    out_png = Path(sys.argv[2])
    title = sys.argv[3]

    if not data_file.exists():
        print("Plik nie istnieje:", data_file)
        return

    nodes, elements, nx, ny = load_sections(data_file)

    node_id = nodes[:, 0].astype(int)
    x = nodes[:, 1]
    y = nodes[:, 2]

    plt.figure(figsize=(8, 8))
    plt.scatter(x, y, color="blue", zorder=3)

    # opisy węzłów
    for nid, xn, yn in zip(node_id, x, y):
        plt.annotate(str(nid+1), (xn, yn), textcoords="offset points",
                     xytext=(0, 5), ha="center", color="blue", fontsize=10)

    # rysowanie elementów
    for elem_id, i, j, k in elements:
        # współrzędne węzłów trójkąta
        xs = [x[i], x[j], x[k], x[i]]  # zamknięcie trójkąta
        ys = [y[i], y[j], y[k], y[i]]
        plt.plot(xs, ys, color="gray", zorder=2)

        # środek trójkąta
        cx = (x[i] + x[j] + x[k]) / 3
        cy = (y[i] + y[j] + y[k]) / 3

        # numer elementu w środku
        plt.annotate(str(elem_id+1), (cx, cy), ha="center", va="center",
                     fontsize=12, color="red")

    # ustawienie etykiet osi co 1
    plt.xticks(np.arange(np.floor(min(x)), np.ceil(max(x))+1, 1))
    plt.yticks(np.arange(np.floor(min(y)), np.ceil(max(y))+1, 1))

    plt.title(title)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

    print("Zapisano wykres:", out_png)

if __name__ == "__main__":
    main()