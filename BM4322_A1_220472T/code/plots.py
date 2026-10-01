# =============================================================================
#  plots.py - all the figures of my report, one function per figure.
#  Each function returns the figure, so I can show it in the notebook and save
#  it with save_all_figures() (the report uses these exact file names).
# =============================================================================
from collections import Counter

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties
import matplotlib.transforms as transforms

from config import UPSTREAM, WINDOW, BASES, WEAK, MAIN_T, BLUE, ORANGE, GREEN, GREY
from helpers import result_file, upstream_position

X_POS = np.arange(-UPSTREAM, 0)                 # positions -100 ... -1


def new_figure(width, height, ncols=1, nrows=1, **kwargs):
    plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans"})
    return plt.subplots(nrows, ncols, figsize=(width, height), dpi=200, facecolor="white", **kwargs)


def nice_axes(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.yaxis.grid(True, color="#e1e0d9", linewidth=0.8)
    ax.set_axisbelow(True)


def q1_positions(q1, q5=None):
    """Q1: where the first WWWW was found (with the shuffled curve if q5 is given)."""
    fig, ax = new_figure(7, 3)
    nice_axes(ax)
    ax.bar(X_POS, [q1["q1_counts"].get(-x, 0) for x in X_POS], width=0.8, color=BLUE, label="Upstream sequences")
    if q5 is not None:
        shuf = q5["control"]["q1_pos"]
        shuf_counts = Counter(shuf)
        ax.plot(X_POS, np.array([shuf_counts.get(-x, 0) for x in X_POS]) / len(shuf) * len(q1["found"]),
                color=GREY, lw=1.4, label="Shuffled sequences (chance, scaled)")
    ax.set_xlim(-101, 0)
    ax.set_xlabel("Position of the WWWW hit relative to the start codon (bases)")
    ax.set_ylabel("Number of genes")
    ax.set_title(f"Q1 - first intact WWWW found by the local search (n = {len(q1['found'])})", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    return fig


def all_occurrences(d, q1):
    """All WWWW occurrences (not only the first) and the A/T content at each position."""
    occ_counts = q1["occ_counts"]
    at_profile = np.array([[b in WEAK for b in seq] for seq in d["upstream"]]).mean(axis=0) * 100
    p_w, gc = d["p_w"], d["gc"]
    fig, (a1, a2) = new_figure(7, 4.4, nrows=2, sharex=True)
    nice_axes(a1)
    nice_axes(a2)
    a1.bar(X_POS, [occ_counts.get(-x, 0) / 10 for x in X_POS], width=0.8, color=GREEN)
    a1.set_ylabel("% of sequences\nwith WWWW here")
    a1.set_title("All WWWW occurrences (about 12 per sequence) and the A/T content at each position",
                 loc="left", fontsize=10)
    a2.plot(X_POS, at_profile, color=ORANGE, lw=1.6)
    a2.axhline(100 * p_w, color=GREY, lw=1)
    a2.axhline(100 * (1 - gc), color="#c3c2b7", lw=1)
    box = dict(facecolor="white", edgecolor="none", pad=1)
    a2.text(-58, 64, f"mean of the upstream regions = {100 * p_w:.1f}% (grey line)", fontsize=8, color=GREY, bbox=box)
    a2.text(-99, 100 * (1 - gc) - 5.5, f"whole genome = {100 * (1 - gc):.1f}%", fontsize=8, color=GREY, bbox=box)
    a2.set_ylim(30, 68)
    a2.set_ylabel("A/T at this position (%)")
    a2.set_xlabel("Position relative to the start codon (bases)")
    a2.set_xlim(-101, 0)
    fig.tight_layout()
    return fig


def q2_runs(d, q2):
    """Q2: number of consecutive Ws, observed and expected for random DNA."""
    run_counts, expected, p_w = q2["run_counts"], q2["expected"], d["p_w"]
    fig, ax = new_figure(6, 3)
    nice_axes(ax)
    labels = list(range(4, 15)) + ["≥15"]
    obs = [run_counts.get(n, 0) for n in range(4, 15)] + [q2["obs_ge15"]]
    exp = [expected[n] for n in range(4, 15)] + [q2["exp_ge15"]]
    ax.bar(range(len(labels)), obs, width=0.6, color=BLUE, label="Observed")
    ax.plot(range(len(labels)), exp, color=ORANGE, marker="o", ms=4, lw=1.6,
            label=f"Expected for random bases (P(W) = {p_w:.3f})")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels([str(v) for v in labels])
    ax.set_xlabel("Number of consecutive Ws at the hit")
    ax.set_ylabel("Number of genes")
    ax.text(0.98, 0.5, f"8 or more Ws:\nobserved {q2['obs_ge8']}, expected {q2['exp_ge8']:.0f}", transform=ax.transAxes,
            ha="right", va="top", fontsize=8, color=GREY)
    ax.set_title(f"Q2 - length of the A/T run at each hit (n = {len(q2['runs'])})", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    return fig


def draw_logo(ax, ppm, info):
    """Sequence logo: the height of each letter = probability x information (bits)."""
    colours = {"A": "#008300", "C": "#2a78d6", "G": "#eda100", "T": "#e34948"}
    font = FontProperties(family="DejaVu Sans", weight="bold")
    for j in range(ppm.shape[1]):
        y = 0.0
        for b in np.argsort(ppm[:, j] * info[j]):          # smallest letter at the bottom
            h = ppm[b, j] * info[j]
            if h < 1e-3:
                continue
            letter = TextPath((0, 0), BASES[b], size=1, prop=font)
            box_ = letter.get_extents()
            t = (transforms.Affine2D().translate(-box_.x0, -box_.y0)
                 .scale(0.9 / box_.width, h / box_.height).translate(j + 0.55, y))
            ax.add_patch(PathPatch(t.transform_path(letter), color=colours[BASES[b]], lw=0))
            y += h
    ax.set_xlim(0.5, ppm.shape[1] + 0.5)
    ax.set_ylim(0, 2)
    ax.set_xticks(range(1, ppm.shape[1] + 1))
    ax.set_xlabel("Position in the 6-base window")
    ax.set_ylabel("Information (bits)")


def q3_ppm_logo(q3):
    """Q3: the PPM as a coloured table and the sequence logo."""
    ppm, info = q3["ppm"], q3["info"]
    fig, (a1, a2) = new_figure(7.2, 2.9, ncols=2, gridspec_kw={"width_ratios": [1.15, 1]})
    a1.imshow(ppm, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    for b in range(4):
        for j in range(WINDOW):
            a1.text(j, b, f"{ppm[b, j]:.3f}", ha="center", va="center", fontsize=8,
                    color="white" if ppm[b, j] > 0.45 else "black")
    a1.set_xticks(range(WINDOW))
    a1.set_xticklabels([str(j + 1) for j in range(WINDOW)])
    a1.set_yticks(range(4))
    a1.set_yticklabels(list(BASES))
    a1.set_xlabel("Position in the 6-base window")
    a1.set_title("Position probability matrix", loc="left", fontsize=10)
    nice_axes(a2)
    a2.yaxis.grid(False)
    draw_logo(a2, ppm, info)
    a2.set_title("Sequence logo", loc="left", fontsize=10)
    fig.tight_layout()
    return fig


def q4_positions(q4, q5=None):
    """Q4: position of the best-scoring window (with the shuffled curve if q5 is given)."""
    det_pos = q4["det_pos"]
    fig, ax = new_figure(7, 3)
    nice_axes(ax)
    det_counts = Counter(det_pos.tolist())
    ax.bar(X_POS, [det_counts.get(-x, 0) for x in X_POS], width=0.8, color=ORANGE, label="Upstream sequences")
    if q5 is not None:
        shuf = q5["control"]["det_pos"]
        shuf_counts = Counter(shuf.tolist())
        xs = X_POS[X_POS <= -WINDOW]
        ax.plot(xs, np.array([shuf_counts.get(-x, 0) for x in xs]) / len(shuf) * len(det_pos),
                color=GREY, lw=1.4, label="Shuffled sequences (chance, scaled)")
    ax.set_xlim(-101, 0)
    ax.set_ylim(0, max(det_counts.values()) * 1.32)
    ax.set_xlabel("Position of the best-scoring 6-base window relative to the start codon (bases)")
    ax.set_ylabel("Number of genes")
    ax.set_title(f"Q4 - statistical alignment, detected at T = {MAIN_T} (n = {len(det_pos)})", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8, loc="upper left", ncol=2)
    fig.tight_layout()
    return fig


def q5_threshold(q1, q5):
    """Q5: detection rate against the threshold, real and shuffled."""
    n_found = len(q1["found"])
    fig, ax = new_figure(6, 3)
    nice_axes(ax)
    ax.plot(q5["sweep"], q5["real_curve"], color=BLUE, lw=2, label="Upstream sequences (statistical alignment)")
    ax.plot(q5["sweep"], q5["shuf_curve"], color=ORANGE, lw=2, label="Shuffled sequences")
    ax.axhline(n_found / 10, color=GREEN, lw=1.6, label=f"Local search ({n_found / 10:.1f}%)")
    ax.set_xlim(0.1, -12.1)
    ax.set_ylim(0, 102)
    ax.set_xlabel("Threshold = best score - consensus score (ln units)")
    ax.set_ylabel("Genes detected (%)")
    ax.set_title("Q4/Q5 - detection rate against the threshold", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout()
    return fig


def q5_positions(q1, q4, q5):
    """Q5: the two position distributions, and the two positions of each gene."""
    fig, (a1, a2) = new_figure(7.4, 3.1, ncols=2)
    nice_axes(a1)
    nice_axes(a2)
    for values, colour, label in ((q1["q1_pos"], BLUE, "Local search (Q1)"),
                                  (q4["det_pos"], ORANGE, f"Statistical alignment (Q4, T = {MAIN_T})")):
        v = np.sort(-np.asarray(values))
        a1.step(v, np.arange(1, len(v) + 1) / len(v) * 100, where="post", color=colour, lw=2, label=label)
    a1.set_xlim(-101, 0)
    a1.set_xlabel("Position (bases)")
    a1.set_ylabel("Cumulative % of genes")
    a1.set_title("Position distributions", loc="left", fontsize=10)
    a1.legend(frameon=False, fontsize=7.5, loc="lower right")
    rng = np.random.default_rng(1)                 # small jitter so that points do not hide each other
    a2.scatter([-upstream_position(q1["hit_start"][k]) + rng.uniform(-0.3, 0.3) for k in q5["both"]],
               [-q4["stat_pos"][k] + rng.uniform(-0.3, 0.3) for k in q5["both"]],
               s=6, color=BLUE, alpha=0.45, linewidths=0)
    a2.plot([-100, 0], [-100, 0], color="#898781", lw=1)
    a2.set_xlim(-101, 0)
    a2.set_ylim(-101, 0)
    a2.set_xlabel("Local search position")
    a2.set_ylabel("Statistical alignment position")
    a2.set_title("Same gene, two methods", loc="left", fontsize=10)
    fig.tight_layout()
    return fig


def save_all_figures(d, q1, q2, q3, q4, q5):
    """Make the 7 report figures, save them in results/figures and close them."""
    figures = {"fig_Q1_positions.png": q1_positions(q1, q5),
               "fig_Q1_all_occurrences.png": all_occurrences(d, q1),
               "fig_Q2_consecutive_W.png": q2_runs(d, q2),
               "fig_Q3_PPM_logo.png": q3_ppm_logo(q3),
               "fig_Q4_positions.png": q4_positions(q4, q5),
               "fig_Q5_threshold.png": q5_threshold(q1, q5),
               "fig_Q5_positions.png": q5_positions(q1, q4, q5)}
    for name, fig in figures.items():
        fig.savefig(result_file("figures", name), facecolor="white")
        plt.close(fig)
    return list(figures)
