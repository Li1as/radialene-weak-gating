"""Regenerate radialene_figure.pdf with matplotlib and numpy."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def tuning_curve(q, y):
    older = np.ones_like(y)
    latest = np.ones_like(y)
    for _ in range(2, q + 1):
        older, latest = latest, latest - y*y*older
    return y**q / latest


def main():
    plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans",
                         "pdf.fonttype": 42, "axes.spines.top": False,
                         "axes.spines.right": False})
    fig, (skeleton, curves) = plt.subplots(1, 2, figsize=(6.6, 2.65),
                                           gridspec_kw={"width_ratios": [1, 1.4]})
    colors = {"x": "#a34422", "y": "#176d75", "contact": "#243447"}
    theta = np.pi/2 + (np.arange(5)-1)*2*np.pi/5
    ring = np.stack((np.cos(theta), np.sin(theta)), axis=1)
    for i, point in enumerate(ring):
        j = (i+1) % 5
        skeleton.plot(*np.stack((point, ring[j])).T, color="0.32", lw=1.5)
        perp = np.array([-point[1], point[0]])*.045
        for offset in (perp, -perp):
            segment = np.stack((point*1.02+offset, point*1.55+offset))
            skeleton.plot(*segment.T, color="0.32", lw=1)
        key = ("contact", "x", "contact", "y", "y")[i]
        label = (r"$s$", r"$x$", r"$t$", r"$y$", r"$y$")[i]
        skeleton.scatter(*(point*1.65), color=colors[key], s=34, zorder=3)
        skeleton.text(*(point*1.99), label, va="center", ha="center",
                      color=colors[key], fontsize=12)
    skeleton.set(xlim=(-2.15,2.15), ylim=(-2.15,2.15), aspect="equal")
    skeleton.axis("off")
    skeleton.set_title("(a)  [5]radialene", loc="left", fontsize=10, pad=2)
    y = np.geomspace(.014, .32, 201)
    for q, color, label in [(2, colors["x"], "[5]radialene"),
                            (4, colors["y"], "[7]radialene")]:
        curves.loglog(y, tuning_curve(q, y), color=color, lw=2, label=label)
        curves.scatter(.1, tuning_curve(q, np.array([.1]))[0], s=18, color=color, zorder=3)
    curves.loglog(y, y, color="0.45", ls="--", lw=1.15, label=r"equal shifts $x=y$")
    curves.set(xlabel=r"Long-arc shift $y$", ylabel=r"Short-arc shift $x$",
               xlim=(.012, .36), ylim=(2e-8, .5))
    curves.set_title("(b)  Exact zero-transmission curves", loc="left", fontsize=10, pad=5)
    curves.legend(frameon=False, loc="lower right", fontsize=8, handlelength=1.6)
    curves.grid(which="major", alpha=.18, lw=.6)
    fig.subplots_adjust(left=.03, right=.98, bottom=.2, top=.89, wspace=.31)
    fig.savefig(Path(__file__).with_name("radialene_figure.pdf"), bbox_inches="tight")


if __name__ == "__main__":
    main()
