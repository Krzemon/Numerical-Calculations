#!/usr/bin/env python3
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

def main():
    if len(sys.argv)!=4:
        print("Uzycie: python3 plot_zad1.py <plik_danych> <plik_wykresu> <tytul>")
        return

    data_file=Path(sys.argv[1])
    out_png=Path(sys.argv[2])
    title=sys.argv[3]

    if not data_file.exists():
        print("Plik nie istnieje:", data_file)
        return

    data = np.loadtxt(data_file)
    plt.figure(figsize=(8,6))
    plt.imshow(data, origin='upper', cmap='viridis', aspect='auto')
    plt.colorbar(label='Wartość')
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    print("Zapisano wykres:", out_png)

if __name__=="__main__":
    main()