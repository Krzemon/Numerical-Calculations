#!/usr/bin/env python3
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

def main():
    if len(sys.argv) != 3:
        print("Użycie: python3 plot_chart.py <plik_danych> <plik_wykresu>")
        return

    data_file = Path(sys.argv[1])
    out_png = Path(sys.argv[2])

    if not data_file.exists():
        print("Plik nie istnieje:", data_file)
        return

    # wczytanie danych
    data = np.loadtxt(data_file)
    t = data[:,0]
    yk_O_c2 = data[:,1]
    yk_O_c3 = data[:,2]
    yk_O_yk = data[:,3]
    yk_E_yk = data[:,4]

    fig, axs = plt.subplots(2, 2, figsize=(10,8))
    
    axs[0,0].plot(t, yk_O_c2, marker='o', markersize=3, linestyle='-', color='blue')
    axs[0,0].set_title("y^T O c2")
    axs[0,0].set_xlabel("t")
    axs[0,0].set_ylabel("y^T O c2")
    axs[0,0].grid(True)

    axs[0,1].plot(t, yk_O_c3, marker='o', markersize=3, linestyle='-', color='green')
    axs[0,1].set_title("y^T O c3")
    axs[0,1].set_xlabel("t")
    axs[0,1].set_ylabel("y^T O c3")
    axs[0,1].grid(True)

    axs[1,0].plot(t, yk_O_yk, marker='o', markersize=3, linestyle='-', color='red')
    axs[1,0].set_title("y^T O y")
    axs[1,0].set_xlabel("t")
    axs[1,0].set_ylabel("y^T O y")
    axs[1,0].grid(True)

    axs[1,1].plot(t, yk_E_yk, marker='o', markersize=3, linestyle='-', color='purple')
    axs[1,1].set_title("y^T E y")
    axs[1,1].set_xlabel("t")
    axs[1,1].set_ylabel("y^T E y")
    axs[1,1].grid(True)

    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    print("Zapisano wykres:", out_png)

if __name__ == "__main__":
    main()