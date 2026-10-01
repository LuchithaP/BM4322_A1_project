# =============================================================================
#  q1_local_search.py - QUESTION 1.
#  Standard local search (Smith-Waterman) for the intact query WWWW in the
#  100 upstream bases of every gene. I also count all WWWW occurrences and look
#  at the genes where no WWWW was found.
# =============================================================================
import re
from collections import Counter

import numpy as np
import pandas as pd

from config import QUERY, MATCH, MISMATCH, GAP, WEAK, UPSTREAM
from helpers import report, heading, upstream_position, result_file


def score(query_letter, base):
    """+1 if the base fits the query letter (W fits A and T), otherwise -1."""
    allowed = WEAK if query_letter == "W" else query_letter
    return MATCH if base in allowed else MISMATCH


def local_search(seq, query):
    """Smith-Waterman local alignment of the query against seq.
    Returns (best score, index in seq where the best alignment starts).
    The query is 'intact' when the best score is len(query) * MATCH, which can
    only happen when every query letter is matched with no gap."""
    n, m = len(seq), len(query)
    H = np.zeros((n + 1, m + 1))
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            H[i, j] = max(0,
                          H[i - 1, j - 1] + score(query[j - 1], seq[i - 1]),   # match / mismatch
                          H[i - 1, j] + GAP,                                   # gap in the query
                          H[i, j - 1] + GAP)                                   # gap in the sequence
    best = H.max()
    # np.argmax gives the FIRST cell with the best score, so when there are many
    # equal hits I get the one nearest to the 5' end (the same as max() in MATLAB).
    i, j = np.unravel_index(np.argmax(H), H.shape)
    while i > 0 and j > 0 and H[i, j] > 0:        # trace back to the start of the alignment
        if H[i, j] == H[i - 1, j - 1] + score(query[j - 1], seq[i - 1]):
            i, j = i - 1, j - 1
        elif H[i, j] == H[i - 1, j] + GAP:
            i -= 1
        else:
            j -= 1
    return best, i


def search_all(upstream):
    """Run the local search on every sequence. Returns scores, hit starts and intact flags."""
    sw_score, hit_start = [], []
    for seq in upstream:
        s, start = local_search(seq, QUERY)
        sw_score.append(s)
        hit_start.append(int(start))
    intact = [s == len(QUERY) * MATCH for s in sw_score]
    # check: an intact hit must be the first run of 4 A/T bases (regular expression)
    for seq, ok, start in zip(upstream, intact, hit_start):
        m = re.search("[AT]{4}", seq)
        assert ok == (m is not None)
        if ok:
            assert start == m.start()
    return sw_score, hit_start, intact


def genes_without_hit(d, sw_score, intact):
    """For the genes with no intact WWWW: G+C content and how much of the upstream
    region is inside another gene."""
    genes, upstream = d["genes"], d["upstream"]
    rows = []
    for k in range(len(genes)):
        if not intact[k]:
            b = genes.loc[k, "Begin"]
            gc_up = (upstream[k].count("G") + upstream[k].count("C"))
            inside = d["coverage"][b - UPSTREAM:b].mean()
            rows.append([genes.loc[k, "Locus tag"], genes.loc[k, "Symbol"], sw_score[k], gc_up,
                         round(100 * inside), upstream[k]])
            report(f"  no hit: {genes.loc[k, 'Locus tag']} {genes.loc[k, 'Symbol']}  best score {sw_score[k]:.0f}, "
                   f"GC {gc_up}%, {100 * inside:.0f}% of the upstream bases are inside another gene")
    pd.DataFrame(rows, columns=["locus_tag", "gene", "best_score", "GC_percent",
                                "percent_inside_other_gene", "upstream_100nt"]).to_csv(
        result_file("tables", "Q1_genes_without_WWWW.csv"), index=False)
    return rows


def run_q1(d):
    """Question 1. d is the dictionary from run_step0()."""
    heading("Q1 - standard local search for WWWW")
    genes, upstream = d["genes"], d["upstream"]
    sw_score, hit_start, intact = search_all(upstream)
    found = [k for k in range(len(genes)) if intact[k]]
    q1_pos = [upstream_position(hit_start[k]) for k in found]
    report(f"Best score distribution: { {int(k): v for k, v in Counter(sw_score).items()} }")
    report(f"Genes with an intact WWWW: {len(found)}/{len(genes)} ({100 * len(found) / len(genes):.1f}%)")
    report(f"Upstream position of the hit: median -{np.median(q1_pos):.0f}, "
           f"IQR -{np.percentile(q1_pos, 75):.0f} to -{np.percentile(q1_pos, 25):.0f}, "
           f"{q1_pos.count(100)} hits at -100, {sum(u >= 80 for u in q1_pos)} hits at -80 or further")

    # all occurrences of WWWW (to understand the result, not only the first one)
    all_hits = [[m.start() for m in re.finditer("(?=[AT]{4})", seq)] for seq in upstream]
    n_occ = [len(h) for h in all_hits]
    occ_counts = Counter(UPSTREAM - k for h in all_hits for k in h)     # how many hits start at each position
    last_pos = [upstream_position(h[-1]) for h in all_hits if h]
    report(f"WWWW occurrences per sequence: mean {np.mean(n_occ):.1f}, median {np.median(n_occ):.0f}, "
           f"range {min(n_occ)}-{max(n_occ)}")
    report(f"If I took the LAST hit (closest to the gene) instead: median -{np.median(last_pos):.0f}")
    no_hit_rows = genes_without_hit(d, sw_score, intact)

    # tables for Q1
    pd.DataFrame({"locus_tag": genes["Locus tag"], "gene": genes["Symbol"], "sw_score": sw_score,
                  "intact_hit": [int(x) for x in intact],
                  "upstream_position": [upstream_position(hit_start[k]) if intact[k] else "" for k in range(len(genes))],
                  "matched_bases": [upstream[k][hit_start[k]:hit_start[k] + 4] if intact[k] else "" for k in range(len(genes))],
                  "n_WWWW_in_sequence": n_occ}).to_csv(result_file("tables", "Q1_local_search_results.csv"), index=False)
    q1_counts = Counter(q1_pos)
    pd.DataFrame({"position": [-u for u in range(UPSTREAM, 0, -1)],
                  "n_genes": [q1_counts.get(u, 0) for u in range(UPSTREAM, 0, -1)]}).to_csv(
        result_file("tables", "Q1_position_distribution.csv"), index=False)

    return {"sw_score": sw_score, "hit_start": hit_start, "intact": intact, "found": found,
            "q1_pos": q1_pos, "q1_counts": q1_counts, "all_hits": all_hits, "n_occ": n_occ, "occ_counts": occ_counts,
            "last_pos": last_pos, "no_hit_rows": no_hit_rows}
