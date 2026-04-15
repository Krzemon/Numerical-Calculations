import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

# ===================== STYL =====================
rc("text.latex", preamble=r"\usepackage{lmodern} \usepackage{physics}")

plt.rcParams.update({
    "font.size": 16,
    "axes.labelsize": 20,
    "figure.figsize": (7, 6),
    "text.usetex": True,
})

# ===================== WCZYTYWANIE =====================
def load_density(filename):
    units = "nm"
    rows = []

    with open(filename) as f:
        for line in f:
            s = line.strip()

            if not s:
                continue

            if s.startswith("#"):
                if "units=au" in s:
                    units = "au"
                elif "units=nm" in s:
                    units = "nm"
                continue

            parts = s.split()
            if len(parts) != 3:
                continue

            rows.append([float(parts[0]), float(parts[1]), float(parts[2])])

    data = np.array(rows)

    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]

    ux = np.unique(x)
    uy = np.unique(y)
    ux.sort()
    uy.sort()
    nx = len(ux)
    ny = len(uy)

    X, Y = np.meshgrid(ux, uy, indexing="ij")
    Z = np.zeros((nx, ny))

    idx_x = {v: i for i, v in enumerate(ux)}
    idx_y = {v: i for i, v in enumerate(uy)}

    for xi, yi, zi in zip(x, y, z):
        i = idx_x[xi]
        j = idx_y[yi]
        Z[i, j] = zi

    return X, Y, Z, units


# ===================== INTERPOLACJA =====================
def interpolate_grid(X, Y, Z, n_new):

    x = X[:, 0]
    y = Y[0, :]

    nx = len(x)
    ny = len(y)
    x_new = np.linspace(x.min(), x.max(), n_new)
    y_new = np.linspace(y.min(), y.max(), n_new)

    Xn, Yn = np.meshgrid(x_new, y_new, indexing="ij")

    # interpolacja po x
    Zn = np.zeros((n_new, ny))

    for j in range(ny):
        Zn[:, j] = np.interp(x_new, x, Z[:, j])

    # interpolacja po y
    Zfinal = np.zeros((n_new, n_new))

    for i in range(n_new):
        Zfinal[i, :] = np.interp(y_new, y, Zn[i, :])

    return Xn, Yn, Zfinal


# ===================== OPIS OSI =====================
def labels(units):
    if units == "nm":
        return (r"$x_1\;[\mathrm{nm}]$",
                r"$x_2\;[\mathrm{nm}]$",
                r"$|\Psi|^2\;[\mathrm{nm}^{-2}]$")
    else:
        return (r"$x_1\;[a_0]$",
                r"$x_2\;[a_0]$",
                r"$|\Psi|^2\;[a_0^{-2}]$")


# ===================== WYKRES =====================
def plot_density(X, Y, Z, units,
                 output_file,
                 title,
                 mode="raw",
                 n_interp=200):

    if mode == "interp":
        X, Y, Z = interpolate_grid(X, Y, Z, n_interp)

    fig, ax = plt.subplots()

    im = ax.pcolormesh(
        X,
        Y,
        Z,
        shading="auto",
        cmap="viridis"
    )

    cbar = fig.colorbar(im, ax=ax)

    xlabel, ylabel, zlabel = labels(units)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    cbar.set_label(zlabel)
    ax.set_aspect("equal")
    fig.savefig(output_file)

    plt.close(fig)
    print(f"Zapisano: {output_file}")


# ===================== MAIN =====================
def main():

    if len(sys.argv) < 3:
        print(
            "Uzycie:\n"
            "python3 plot_density.py <plik_danych> <plik_wyjsciowy> "
            "[--mode raw|interp] "
            "[--n 200]"
        )
        return

    inp = Path(sys.argv[1])
    out = Path(sys.argv[2])

    mode = "raw"
    n_interp = 200

    args = sys.argv[3:]

    if "--mode" in args:
        i = args.index("--mode")
        mode = args[i + 1]

    if "--n" in args:
        i = args.index("--n")
        n_interp = int(args[i + 1])

    if not inp.exists():
        print("Brak pliku:", inp)
        return

    X, Y, Z, units = load_density(inp)

    title = r"$|\Psi(x_1,x_2)|^2$"
    plot_density(X, Y, Z, units, out,title, mode=mode, n_interp=n_interp)


if __name__ == "__main__":
    main()