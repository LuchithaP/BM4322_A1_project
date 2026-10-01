# =============================================================================
#  summary.py - collects every number I used in my report into one dictionary,
#  which is saved as results/summary.json. I filled in the report from this file.
# =============================================================================
from collections import Counter

import numpy as np

from config import UPSTREAM, THRESHOLDS


def pct(values, q):
    return float(np.percentile(values, q))


def build_summary(d, met, q1, q2, q3, q4, q5):
    """d, met, q1 ... q5 are the dictionaries returned by the steps."""
    on_chrom, table, genes, gc, p_w = d["on_chrom"], d["table"], d["genes"], d["gc"], d["p_w"]
    q1_pos, n_occ, det_pos, ctrl = q1["q1_pos"], q1["n_occ"], q4["det_pos"], q5["control"]
    inside = q5["inside"]
    return {
        # Step 0
        "genome_length": len(d["genome"]), "gc_percent": 100 * gc, "genes_in_table": len(table),
        "gene_types": {k: int(v) for k, v in table["Gene Type"].value_counts().items()},
        "chromosome_genes": len(on_chrom), "plus_genes": int(sum(on_chrom["Orientation"] == "plus")),
        "minus_genes": int(sum(on_chrom["Orientation"] == "minus")), "plus_protein_coding": len(
            on_chrom[(on_chrom["Orientation"] == "plus") & (on_chrom["Gene Type"] == "protein-coding")]),
        "first_gene_begin": int(genes["Begin"].iloc[0]), "last_gene_end": int(genes["End"].iloc[-1]),
        "start_codons": dict(met["codon_counts"]), "met": met["n_met"], "length_ok": met["n_length_ok"],
        "upstream_AT_percent": 100 * p_w, "genome_AT_percent": 100 * (1 - gc),
        # Q1
        "q1_found": len(q1["found"]), "q1_median": float(np.median(q1_pos)), "q1_q1": pct(q1_pos, 25),
        "q1_q3": pct(q1_pos, 75), "q1_at_100": q1_pos.count(100), "q1_ge80": sum(u >= 80 for u in q1_pos),
        "occ_mean": float(np.mean(n_occ)), "occ_median": float(np.median(n_occ)), "occ_min": min(n_occ),
        "occ_max": max(n_occ), "occ_total": int(sum(n_occ)), "last_hit_median": float(np.median(q1["last_pos"])),
        "occ_percent_by_position": {u: q1["occ_counts"].get(u, 0) / 10 for u in range(UPSTREAM, 3, -1)},
        "no_hit": q1["no_hit_rows"],
        # Q2
        "q2_mean": float(np.mean(q2["runs"])), "q2_median": float(np.median(q2["runs"])),
        "q2_min": min(q2["runs"]), "q2_max": max(q2["runs"]),
        "q2_counts": {int(n): int(v) for n, v in sorted(q2["run_counts"].items())},
        "q2_expected": {int(n): float(v) for n, v in q2["expected"].items()},
        "q2_obs_ge8": q2["obs_ge8"], "q2_exp_ge8": q2["exp_ge8"], "q2_obs_ge15": q2["obs_ge15"], "q2_exp_ge15": q2["exp_ge15"],
        "q2_longest": q2["longest"],
        # Q3
        "q3_n": len(q3["hexamers"]), "q3_distinct": len(set(q3["hexamers"])),
        "q3_top": Counter(q3["hexamers"]).most_common(5),
        "q3_counts": q3["counts"].astype(int).tolist(), "q3_ppm": q3["ppm"].tolist(), "q3_info": q3["info"].tolist(),
        "q3_info_total": float(q3["info"].sum()), "q3_consensus": q3["consensus"],
        "q3_consensus_score": float(q3["consensus_score"]),
        # Q4
        "q4_detected_percent": {str(t): float(q4["detected"][t].mean() * 100) for t in THRESHOLDS},
        "q4_detected_n": {str(t): int(q4["detected"][t].sum()) for t in THRESHOLDS},
        "q4_no_hit_scores": q4["no_hit_scores"],
        "q4_median": float(np.median(det_pos)), "q4_q1": pct(det_pos, 25), "q4_q3": pct(det_pos, 75),
        "q4_top": q4["top_best"], "q4_all_w": q4["all_w"], "q4_worst_wwww": float(q4["worst_wwww"]),
        "q4_cost_C_pos1": float(q3["log_ppm"][1, 0] - q3["log_ppm"][:, 0].max()),
        # Q5
        "q5_both": len(q5["both"]), "q5_same": q5["same"], "q5_close": q5["close"],
        "q5_shift_median": float(np.median(q5["shift"])), "q5_shift_min": int(min(q5["shift"])),
        "q5_ks_D": float(q5["ks"].statistic), "q5_ks_p": float(q5["ks"].pvalue),
        "q5_run_local": float(np.mean(q2["runs"])), "q5_run_stat": float(np.mean(q5["run_at_stat"])),
        "shuf_found": float(ctrl["found"]), "shuf_occ": float(ctrl["occ"]),
        "shuf_det": {str(t): float(v) for t, v in ctrl["det"].items()},
        "shuf_q1_median": float(np.median(ctrl["q1_pos"])), "shuf_q4_median": float(np.median(ctrl["det_pos"])),
        "regions": q5["region_rows"],
        "inside_other_gene": int(inside.sum()), "inside_same_strand_close": q5["same_strand_close"],
        "found_inside_percent": float(q5["found_inside"]), "found_free_percent": float(q5["found_free"]),
        "sweep": {f"{t:.2f}": [float(r), float(s)] for t, r, s in zip(q5["sweep"], q5["real_curve"], q5["shuf_curve"])},
        "expected_hits": q5["expected_hits"], "p_at_least_one": q5["p_at_least_one"],
    }
