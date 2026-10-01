# =============================================================================
#  q5_comparison.py - QUESTION 5.
#  I compare the local search (Q1) with the statistical alignment (Q4): how many
#  genes each one detects and where it puts the promoter. To know what the two
#  searches find just by chance, I also shuffle every sequence 100 times and
#  run both searches again on the shuffled sequences.
# =============================================================================
import re
import random

import numpy as np
import pandas as pd
from scipy import stats

from config import UPSTREAM, N_SHUFFLE, SEED, THRESHOLDS, MAIN_T
from helpers import report, heading, result_file
from q2_consecutive_w import count_consecutive_w
from q4_statistical_alignment import stat_align

REGIONS = {"-100 to -26": (26, 100), "-25 to -14": (14, 25), "-13 to -9": (9, 13), "-8 to -6": (6, 8)}


def compare_positions(d, q1, q4):
    """Same genes, two positions: do the two searches agree?"""
    hit_start, best_index = q1["hit_start"], q4["best_index"]
    both = [k for k in q1["found"] if q4["detected"][MAIN_T][k]]
    same = sum(1 for k in both if hit_start[k] == best_index[k])
    close = sum(1 for k in both if abs(hit_start[k] - best_index[k]) <= 2)
    shift = [best_index[k] - hit_start[k] for k in both]          # >= 0 means closer to the gene
    ks = stats.ks_2samp(q1["q1_pos"], q4["det_pos"])
    run_at_stat = [count_consecutive_w(d["upstream"][k], best_index[k]) for k in both]
    return {"both": both, "same": same, "close": close, "shift": shift, "ks": ks, "run_at_stat": run_at_stat}


def shuffle_sequences(upstream):
    """Every sequence shuffled N_SHUFFLE times: same bases, random order.
    The seed is set here, just before the shuffling, so the result is always the same."""
    random.seed(SEED)
    shuffled = []
    for seq in upstream:
        for _ in range(N_SHUFFLE):
            shuffled.append("".join(random.sample(seq, len(seq))))
    return shuffled


def shuffled_control(d, q3):
    """Run both searches on the shuffled sequences. For 100,000 sequences Smith-Waterman
    in pure Python is too slow, so here I use the regular expression - in Q1 I checked
    that it gives exactly the same hits as the Smith-Waterman search."""
    shuffled = shuffle_sequences(d["upstream"])
    first = [re.search("[AT]{4}", s) for s in shuffled]
    c = {"found": np.mean([m is not None for m in first]) * 100,
         "q1_pos": [UPSTREAM - m.start() for m in first if m is not None],
         "occ": np.mean([len(re.findall("(?=[AT]{4})", s)) for s in shuffled])}
    rel, best = [], []
    for s in shuffled:
        sc = stat_align(s, q3["log_ppm"])
        rel.append(sc.max() - q3["consensus_score"])
        best.append(UPSTREAM - int(np.argmax(sc)))
    c["rel"] = np.array(rel)
    c["best"] = np.array(best)
    c["det"] = {t: (c["rel"] >= t).mean() * 100 for t in THRESHOLDS}
    c["det_pos"] = c["best"][c["rel"] >= MAIN_T]
    return c


def share(pos, lo, hi):
    """Percentage of positions between -hi and -lo."""
    pos = np.asarray(pos)
    return 100 * np.mean((pos >= lo) & (pos <= hi))


def hits_by_region(q1, q4, ctrl):
    rows = []
    for name, (lo, hi) in REGIONS.items():
        rows.append([name, share(q1["q1_pos"], lo, hi), share(ctrl["q1_pos"], lo, hi),
                     share(q4["det_pos"], lo, hi), share(ctrl["det_pos"], lo, hi)])
        report(f"  {name:12s} local {rows[-1][1]:5.1f}% (shuffled {rows[-1][2]:5.1f}%)   "
               f"statistical {rows[-1][3]:5.1f}% (shuffled {rows[-1][4]:5.1f}%)")
    pd.DataFrame(rows, columns=["region", "local_real_%", "local_shuffled_%",
                                "statistical_real_%", "statistical_shuffled_%"]).round(2).to_csv(
        result_file("tables", "Q5_hits_by_region.csv"), index=False)
    return rows


def overlap_with_other_genes(d, q1):
    """Is the upstream region free DNA, or inside the previous gene (like inside an operon)?"""
    genes, on_chrom, intact = d["genes"], d["on_chrom"], q1["intact"]
    inside = np.array([d["coverage"][b - UPSTREAM:b].any() for b in genes["Begin"]])
    same_strand_close = 0
    for k in np.where(inside)[0]:
        b = genes.loc[k, "Begin"]
        before = on_chrom[(on_chrom["Begin"] < b) & (on_chrom["Locus tag"] != genes.loc[k, "Locus tag"])]
        nearest = before.loc[before["End"].idxmax()]          # the gene that ends closest to our gene
        if nearest["Orientation"] == "plus" and b - nearest["End"] - 1 < UPSTREAM:
            same_strand_close += 1
    found_inside = np.mean([intact[k] for k in np.where(inside)[0]]) * 100
    found_free = np.mean([intact[k] for k in np.where(~inside)[0]]) * 100
    report(f"Upstream region overlaps another gene: {inside.sum()} genes "
           f"({same_strand_close} of them right after another + strand gene, like inside an operon)")
    report(f"  WWWW found in {found_inside:.1f}% of these and in {found_free:.1f}% of the free (intergenic) ones")
    return {"inside": inside, "same_strand_close": same_strand_close,
            "found_inside": found_inside, "found_free": found_free}


def threshold_sweep(q4, ctrl):
    """Detection rate for thresholds from 0 to -12 (for the figure)."""
    sweep = np.round(np.arange(0, -12.01, -0.25), 2)
    real_curve = [(q4["best_rel"] >= t).mean() * 100 for t in sweep]
    shuf_curve = [(ctrl["rel"] >= t).mean() * 100 for t in sweep]
    pd.DataFrame({"threshold": sweep, "real_%": np.round(real_curve, 2),
                  "shuffled_%": np.round(shuf_curve, 2)}).to_csv(
        result_file("tables", "Q4_Q5_detection_vs_threshold.csv"), index=False)
    return sweep, real_curve, shuf_curve


def chance_expectation(p_w):
    """Expected number of WWWW in 97 positions of random DNA, and P(at least one)."""
    expected_hits = (UPSTREAM - 3) * p_w ** 4
    state = [1.0, 0.0, 0.0, 0.0]                   # probability of the current A/T run length 0..3
    for _ in range(UPSTREAM):
        new = [0.0, 0.0, 0.0, 0.0]
        for r in range(4):
            new[0] += state[r] * (1 - p_w)         # a G or C breaks the run
            if r + 1 < 4:
                new[r + 1] += state[r] * p_w        # an A or T makes the run longer
        state = new                                # runs reaching 4 are removed (= found)
    return expected_hits, 1 - sum(state)


def run_q5(d, q1, q2, q3, q4):
    """Question 5. Uses the results of all the other steps."""
    heading("Q5 - comparison")
    cmp = compare_positions(d, q1, q4)
    n_both = len(cmp["both"])
    report(f"Detected by both: {n_both}; same position {cmp['same']} ({100 * cmp['same'] / n_both:.1f}%), "
           f"within 2 bases {cmp['close']} ({100 * cmp['close'] / n_both:.1f}%)")
    report(f"Statistical hit is never upstream of the local hit: {min(cmp['shift']) >= 0}; "
           f"median shift towards the gene {np.median(cmp['shift']):.0f} bases")
    report(f"KS test between the two position distributions: D = {cmp['ks'].statistic:.3f}, "
           f"p = {cmp['ks'].pvalue:.2e}")
    report(f"Mean A/T run at the hit: local {np.mean(q2['runs']):.2f}, statistical {np.mean(cmp['run_at_stat']):.2f}")

    ctrl = shuffled_control(d, q3)
    report(f"Shuffled: intact WWWW in {ctrl['found']:.1f}% (real {100 * len(q1['found']) / len(d['genes']):.1f}%), "
           f"{ctrl['occ']:.1f} WWWW per sequence (real {np.mean(q1['n_occ']):.1f})")
    report(f"Shuffled: statistical detection {', '.join(f'{t}: {v:.1f}%' for t, v in ctrl['det'].items())}")
    report(f"Shuffled: local-search position median -{np.median(ctrl['q1_pos']):.0f}, "
           f"statistical position median -{np.median(ctrl['det_pos']):.0f}")

    region_rows = hits_by_region(q1, q4, ctrl)
    overlap = overlap_with_other_genes(d, q1)
    sweep, real_curve, shuf_curve = threshold_sweep(q4, ctrl)
    expected_hits, p_at_least_one = chance_expectation(d["p_w"])
    report(f"Random bases with P(W) = {d['p_w']:.3f}: expected {expected_hits:.1f} WWWW per sequence, "
           f"P(at least one) = {p_at_least_one:.3f}")

    return {**cmp, "control": ctrl, "region_rows": region_rows, **overlap, "sweep": sweep,
            "real_curve": real_curve, "shuf_curve": shuf_curve,
            "expected_hits": expected_hits, "p_at_least_one": p_at_least_one}
