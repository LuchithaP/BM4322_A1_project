# =============================================================================
#  q4_statistical_alignment.py - QUESTION 4.
#  Statistical alignment of the same 1000 sequences with the PPM from Q3:
#  I slide a 6-base window along each upstream sequence, score every window with
#  the sum of ln(probabilities), keep the best window and call a promoter when
#  the best score is within T of the consensus score.
# =============================================================================
from collections import Counter

import numpy as np
import pandas as pd

from config import WINDOW, BASES, WEAK, THRESHOLDS, MAIN_T, UPSTREAM
from helpers import report, heading, write_fasta, result_file


def to_numbers(seq):
    """A, C, G, T -> 0, 1, 2, 3 (the row numbers of the PPM)."""
    return np.array([BASES.index(b) for b in seq])


def stat_align(seq, log_ppm):
    """Score every 6-base window of seq: S = ln p(base1, pos1) + ... + ln p(base6, pos6).
    I do it for all windows at once with numpy (the same as two for-loops)."""
    x = to_numbers(seq)
    n_windows = len(x) - WINDOW + 1
    scores = np.zeros(n_windows)
    for j in range(WINDOW):
        scores += log_ppm[x[j:j + n_windows], j]
    return scores


def worst_wwww_score(log_ppm):
    """The lowest score (relative to the consensus) that a window starting with WWWW
    can get: A or T at positions 1-4, anything at positions 5-6."""
    return sum(min(log_ppm[[0, 3], j]) - log_ppm[:, j].max() for j in range(4)) + \
        sum(log_ppm[:, j].min() - log_ppm[:, j].max() for j in range(4, WINDOW))


def run_q4(d, q1, q3):
    """Question 4. d is from run_step0(), q1 from run_q1(), q3 from run_q3()."""
    heading("Q4 - statistical alignment")
    genes, upstream = d["genes"], d["upstream"]
    log_ppm, consensus_score = q3["log_ppm"], q3["consensus_score"]
    n = len(genes)

    best_index, best_rel = [], []
    for seq in upstream:
        s = stat_align(seq, log_ppm)          # 95 windows, from -100 to -6
        best_index.append(int(np.argmax(s)))   # first best window (5' end side) if there is a tie
        best_rel.append(s.max() - consensus_score)
    best_index = np.array(best_index)
    best_rel = np.array(best_rel)
    stat_pos = UPSTREAM - best_index
    best_hex = [upstream[k][best_index[k]:best_index[k] + WINDOW] for k in range(n)]

    detected = {t: best_rel >= t for t in THRESHOLDS}
    for t in THRESHOLDS:
        report(f"Threshold {t}: {detected[t].sum()}/{n} detected ({detected[t].mean() * 100:.1f}%)")
    report(f"Best score minus consensus: median {np.median(best_rel):.3f}, min {best_rel.min():.2f}")
    no_hit = [k for k in range(n) if not q1["intact"][k]]
    report(f"Scores of the {len(no_hit)} genes without WWWW: {np.round(sorted(best_rel[no_hit]), 2)}")
    det_pos = stat_pos[detected[MAIN_T]]
    report(f"Position of the best window (T = {MAIN_T}): median -{np.median(det_pos):.0f}, "
           f"IQR -{np.percentile(det_pos, 75):.0f} to -{np.percentile(det_pos, 25):.0f}")
    top_best = Counter(best_hex[k] for k in range(n) if detected[MAIN_T][k]).most_common(5)
    all_w = sum(1 for k in range(n) if detected[MAIN_T][k] and set(best_hex[k]) <= set(WEAK))
    report(f"Most common best windows: {top_best}; windows made only of A/T: {all_w}")
    worst_wwww = worst_wwww_score(log_ppm)
    report(f"Any window starting with WWWW scores at least {worst_wwww:.2f} below the consensus")

    # files for Q4
    table_q4 = pd.DataFrame({"locus_tag": genes["Locus tag"], "gene": genes["Symbol"],
                             "best_position": -stat_pos, "best_6mer": best_hex,
                             "score_minus_consensus": np.round(best_rel, 4)})
    for t in THRESHOLDS:
        table_q4[f"detected_T{t}"] = detected[t].astype(int)
    table_q4.to_csv(result_file("tables", "Q4_statistical_alignment_results.csv"), index=False)
    write_fasta(result_file("sequences", "Q4_best_6mers.fasta"),
                [(f"{genes.loc[k, 'Locus tag']}|pos=-{stat_pos[k]}|score-consensus={best_rel[k]:.3f}", best_hex[k])
                 for k in range(n)])

    return {"best_index": best_index, "best_rel": best_rel, "stat_pos": stat_pos, "best_hex": best_hex,
            "detected": detected, "det_pos": det_pos, "top_best": top_best, "all_w": all_w,
            "worst_wwww": worst_wwww, "no_hit_scores": sorted(best_rel[no_hit].tolist())}
