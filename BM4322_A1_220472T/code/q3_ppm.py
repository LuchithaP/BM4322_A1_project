# =============================================================================
#  q3_ppm.py - QUESTION 3.
#  I take the 6 bases starting at every WWWW hit (whatever their C/G content),
#  count each base at each position and turn the counts into the position
#  probability matrix (PPM). I also work out the consensus and the information
#  content of every position.
# =============================================================================
from collections import Counter

import numpy as np
import pandas as pd

from config import WINDOW, PSEUDO, BASES
from helpers import report, heading, upstream_position, write_fasta, write_lines, result_file


def get_hexamers(sequences, hit_start, found):
    """The 6-base sequence that starts at each hit."""
    return [sequences[k][hit_start[k]:hit_start[k] + WINDOW] for k in found]


def make_ppm(hexamers):
    """Position frequency matrix (counts) and position probability matrix.
    Rows are A, C, G, T and columns are the 6 positions."""
    counts = np.zeros((4, WINDOW))
    for h in hexamers:
        for j, base in enumerate(h):
            counts[BASES.index(base), j] += 1
    ppm = (counts + PSEUDO) / (len(hexamers) + 4 * PSEUDO)     # each column adds up to 1
    return counts, ppm


def information_content(ppm):
    """Information of each column in bits, compared with a random base (1/4 each)."""
    return (ppm * np.log2(ppm / 0.25)).sum(axis=0)


def run_q3(d, q1):
    """Question 3. d is from run_step0(), q1 from run_q1()."""
    heading("Q3 - position probability matrix")
    genes = d["genes"]
    hit_start, found = q1["hit_start"], q1["found"]

    hexamers = get_hexamers(d["sequences"], hit_start, found)
    counts, ppm = make_ppm(hexamers)
    log_ppm = np.log(ppm)
    consensus = "".join(BASES[b] for b in ppm.argmax(axis=0))
    consensus_score = log_ppm.max(axis=0).sum()
    info = information_content(ppm)

    report(f"Number of 6-base sequences: {len(hexamers)} ({len(set(hexamers))} different)")
    report(f"Most common: {Counter(hexamers).most_common(5)}")
    report("Frequency matrix (rows A C G T):\n" + str(counts.astype(int)))
    report("PPM:\n" + str(np.round(ppm, 3)))
    report(f"Consensus {consensus}, consensus score (sum of ln p) = {consensus_score:.3f}")
    report(f"Information content per position (bits): {np.round(info, 3)}  total {info.sum():.2f}")

    # files for Q3
    write_fasta(result_file("sequences", "Q3_hexamers.fasta"),
                [(f"{genes.loc[k, 'Locus tag']}|pos=-{upstream_position(hit_start[k])}", h)
                 for k, h in zip(found, hexamers)])
    write_lines(result_file("sequences", "Q3_hexamers.txt"), hexamers)
    cols = [f"pos{j + 1}" for j in range(WINDOW)]
    pd.DataFrame(counts.astype(int), index=list(BASES), columns=cols).to_csv(
        result_file("tables", "Q3_frequency_matrix.csv"))
    ppm_table = pd.DataFrame(np.round(ppm, 6), index=list(BASES), columns=cols)
    ppm_table.loc["info_bits"] = np.round(info, 4)
    ppm_table.to_csv(result_file("tables", "Q3_PPM.csv"))

    return {"hexamers": hexamers, "counts": counts, "ppm": ppm, "log_ppm": log_ppm,
            "consensus": consensus, "consensus_score": consensus_score, "info": info}
