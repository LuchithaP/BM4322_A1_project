# =============================================================================
#  q2_consecutive_w.py - QUESTION 2.
#  For every WWWW hit I count how many A/T bases follow each other from the start
#  of the hit, and compare the distribution with what random DNA would give.
# =============================================================================
from collections import Counter

import numpy as np
import pandas as pd

from config import WEAK
from helpers import report, heading, upstream_position, result_file


def count_consecutive_w(seq, start):
    """Number of A/T bases in a row, starting at position 'start' of seq."""
    k = start
    while k < len(seq) and seq[k] in WEAK:
        k += 1
    return k - start


def run_q2(d, q1):
    """Question 2. d is from run_step0(), q1 from run_q1()."""
    heading("Q2 - consecutive Ws")
    genes, upstream, p_w = d["genes"], d["upstream"], d["p_w"]
    hit_start, found = q1["hit_start"], q1["found"]

    runs = [count_consecutive_w(upstream[k], hit_start[k]) for k in found]
    run_counts = Counter(runs)
    report(f"Run length: mean {np.mean(runs):.2f}, median {np.median(runs):.0f}, range {min(runs)}-{max(runs)}")
    report(f"Distribution: { {int(n): v for n, v in sorted(run_counts.items())} }")

    # If the bases were random with P(W) = p_w, the run would continue with probability p_w
    # each time, so the length after the first 4 Ws is geometric: P(n) = p_w^(n-4) * (1-p_w)
    expected = {n: len(runs) * p_w ** (n - 4) * (1 - p_w) for n in range(4, max(runs) + 1)}
    obs_ge8 = sum(v for n, v in run_counts.items() if n >= 8)
    exp_ge8 = len(runs) * p_w ** 4
    obs_ge15 = sum(v for n, v in run_counts.items() if n >= 15)        # the last bar of the figure
    exp_ge15 = len(runs) * p_w ** 11
    report(f"Runs of 8 or more: observed {obs_ge8}, expected for random bases {exp_ge8:.0f}")
    longest = [(genes.loc[k, "Symbol"], -int(upstream_position(hit_start[k])),
                upstream[k][hit_start[k]:hit_start[k] + r]) for k, r in zip(found, runs) if r >= 20]
    report(f"Longest runs: {longest}")

    pd.DataFrame({"n_consecutive_W": list(range(4, max(runs) + 1)),
                  "n_genes": [run_counts.get(n, 0) for n in range(4, max(runs) + 1)],
                  "expected_random": [round(expected[n], 2) for n in range(4, max(runs) + 1)]}).to_csv(
        result_file("tables", "Q2_consecutive_W_distribution.csv"), index=False)

    return {"runs": runs, "run_counts": run_counts, "expected": expected,
            "obs_ge8": obs_ge8, "exp_ge8": exp_ge8, "obs_ge15": obs_ge15, "exp_ge15": exp_ge15,
            "longest": longest}
