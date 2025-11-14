#!/usr/bin/env python3
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

def main():
    data_file = Path(sys.argv[1])
    out_png = Path(sys.argv[2])
    title = sys.argv[3]
    if not data_file.exists(): print("Plik nie istnieje:", data_file); return
    data = np.loadtxt(data_file, comments="#")
    x = data[:,0]
    n_curves = min(5, data.shape[1]-1)
    plt.figure(figsize=(9,6))
    for j in range(1,n_curves+1): plt.plot(x,data[:,j], label=f"u{j-1}")
    plt.xlabel("x"); plt.ylabel("u_mu(x)"); plt.title(title)
    plt.legend(); plt.grid(True); plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    print("Zapisano wykres:", out_png)

if __name__=="__main__": main()
