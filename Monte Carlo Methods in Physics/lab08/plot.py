# basic imports
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import os

# ============ for latex fonts ============
from matplotlib import rc  # , font_manager

rc(
    "text.latex", preamble=r"\usepackage{lmodern}"
)  # this helps use the plots in tex files
plt.rcParams.update({"font.size": 11})
plt.rcParams.update(
    {
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.minor.width": 0.8,
        "ytick.minor.width": 0.8,
        "axes.labelsize": 12,
        "axes.titlesize": 13,
        "axes.linewidth": 0.8,
        "lines.linewidth": 0.8,
        "patch.linewidth": 0.8,
        "legend.fontsize": 10,
        "legend.title_fontsize": 10,
        "legend.fancybox": False,
        "legend.frameon": False,
        "legend.fontsize": 10,
        "legend.handlelength": 1.0,
        # "legend.handletextpad": 0.5,
        "legend.labelspacing": 0.4,
        "figure.titlesize": 12,
        "figure.figsize": [3.54, 3.54],
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "mathtext.fontset": "cm",
        "text.usetex": True,
        "font.family": "Computer Modern Roman",
    }
)
# ==========================================

###################

data_dir = Path("build/data")
plot_dir = Path("plot")

if not Path.exists(plot_dir):
    Path.mkdir(plot_dir)

###################


def load_data(filename):
    data = np.loadtxt(data_dir / filename)
    return (col for col in data.T)


def prep_image_data(x, y, z):
    x = np.unique(x)
    y = np.unique(y)

    X, Y = np.meshgrid(x, y)
    Z = z.reshape(len(x), len(y))

    return (np.transpose(Z), [min(x), max(x), min(x), max(x)])


def set_labels(xlabel="x", ylabel="y", title=""):
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)


def reset_plot_setting():
    plt.close("all")
    # plt.figure(figsize=figsize, dpi=dpi)


def save_plot(filename, bbox_inches="tight"):
    plt.savefig(plot_dir / filename, bbox_inches=bbox_inches)
    reset_plot_setting()


##########################

# rysowanie fullerenow


def plot_fullerene(filename, r_max=2.0, axes_extend=4):

    x, y, z = np.loadtxt(data_dir / filename, unpack=True)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # Plot vertices
    ax.scatter(x, y, z, c="r", marker="o")

    # Plot edges
    for i in range(len(x)):
        for j in range(i + 1, len(x)):
            dist = np.sqrt((x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2 + (z[i] - z[j]) ** 2)
            if dist <= r_max:
                ax.plot([x[i], x[j]], [y[i], y[j]], [z[i], z[j]], c="b")

    ax.set_xlim(-axes_extend, axes_extend)
    ax.set_ylim(-axes_extend, axes_extend)
    ax.set_zlim(-axes_extend, axes_extend)

    ax.set_box_aspect([1, 1, 1])

    return fig, ax


# zad2 - zad5
for i in [2, 3, 4, 5]:
    if i == 2:
        plot_fullerene(f"zad{i}_C60.dat", 0.5 * 3.5244, np.ceil(3.5244))
        save_plot(f"zad{i}_C60.png")
    else:
        it, beta, V_tot, r_avg = load_data(f"zad{i}.dat")

        # E and beta
        fig, ax1 = plt.subplots()

        color = "black"
        ax1.set_xlabel("Iteration number")
        ax1.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
        ax1.set_ylim(-450, -250)
        ax1.set_ylabel(r"$E$", color=color)
        ax1.plot(it, V_tot, color=color)
        ax1.tick_params(axis="y", labelcolor=color)

        ax2 = ax1.twinx()

        color = "red"
        ax2.set_ylabel(r"$\beta$", color=color)
        ax2.plot(it, beta, color=color)
        ax2.tick_params(axis="y", labelcolor=color)

        save_plot(f"zad{i}_E.png")

        # r_avg and beta
        fig, ax1 = plt.subplots()

        color = "black"
        ax1.set_xlabel("Iteration number")
        ax1.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
        ax1.set_ylabel(r"$\overline{r}$", color=color)
        ax1.plot(it, r_avg, color=color)
        ax1.tick_params(axis="y", labelcolor=color)

        ax2 = ax1.twinx()

        color = "red"
        ax2.set_ylabel(r"$\beta$", color=color)
        ax2.plot(it, beta, color=color)
        ax2.tick_params(axis="y", labelcolor=color)

        save_plot(f"zad{i}_r.png")

        # fullerene
        plot_fullerene(f"zad{i}_C60.dat", 0.5 * max(r_avg), np.ceil(max(r_avg)))
        save_plot(f"zad{i}_C60.png")

    # pcf
    pcf = load_data(f"zad{i}_PCF.dat")
    plt.plot(np.linspace(0, 2.5, 100), list(pcf))
    set_labels(
        r"$r$ in units of $\overline{r}$",
        r"$pcf(r)$",
    )
    save_plot(f"zad{i}_PCF.png")


# zad6
for i in range(6):
    it, beta, V_tot, r_avg = load_data(f"zad6_{i}.dat")

    # E and beta
    fig, ax1 = plt.subplots()

    color = "black"
    ax1.set_xlabel("Iteration number")
    ax1.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
    ax1.set_ylim(-450, -250)
    ax1.set_ylabel(r"$E$", color=color)
    ax1.plot(it, V_tot, color=color)
    ax1.tick_params(axis="y", labelcolor=color)

    ax2 = ax1.twinx()

    color = "red"
    ax2.set_ylabel(r"$\beta$", color=color)
    ax2.plot(it, beta, color=color)
    ax2.tick_params(axis="y", labelcolor=color)

    save_plot(f"zad6_E_{i}.png")

    # r_avg and beta
    fig, ax1 = plt.subplots()

    color = "black"
    ax1.set_xlabel("Iteration number")
    ax1.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
    ax1.set_ylabel(r"$\overline{r}$", color=color)
    ax1.plot(it, r_avg, color=color)
    ax1.tick_params(axis="y", labelcolor=color)

    ax2 = ax1.twinx()

    color = "red"
    ax2.set_ylabel(r"$\beta$", color=color)
    ax2.plot(it, beta, color=color)
    ax2.tick_params(axis="y", labelcolor=color)

    save_plot(f"zad6_r_{i}.png")

    # pcf
    pcf = load_data(f"zad6_PCF_{i}.dat")
    plt.plot(np.linspace(0, 2.5, 100), list(pcf))
    set_labels(
        r"$r$ in units of $\overline{r}$",
        r"$pcf(r)$",
    )
    save_plot(f"zad6_PCF_{i}.png")

    # fullerene
    plot_fullerene(f"zad6_C60_{i}.dat", 0.5 * max(r_avg), np.ceil(max(r_avg)))
    save_plot(f"zad6_C60_{i}.png")

# zad7
n, Eb, Eb_std, r_avg, r_avg_std = load_data(f"zad7.dat")

fig, ax = plt.subplots(figsize=[6, 3.54])
ax.plot(n, Eb)
ax.errorbar(n, Eb, yerr=Eb_std, capsize=2, fmt="o", markersize=1.5)
ax.set_xticks([i for i in n if i % 2 == 0])
# add grid x
ax.xaxis.grid(True)
set_labels(r"$n$ - fullerene size", r"$E_b=V_{tot}/n$")
save_plot(f"zad7_E.png")

fig, ax = plt.subplots(figsize=[6, 3.54])
ax.plot(n, r_avg)
ax.errorbar(n, r_avg, yerr=r_avg_std, capsize=2, fmt="o", markersize=1.5)
ax.set_xticks([i for i in n if i % 2 == 0])
# add grid x
ax.xaxis.grid(True)
set_labels(r"$n$ - fullerene size", r"$\overline{r}$")
save_plot(f"zad7_r.png")
