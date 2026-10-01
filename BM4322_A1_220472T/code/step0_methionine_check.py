# =============================================================================
#  step0_methionine_check.py - STEP 0 check of my code.
#  If my positions are right, the 3 downstream bases are a start codon, which
#  codes for Methionine (M). I check a few known genes first, then all 1000.
# =============================================================================
from collections import Counter

import pandas as pd

from config import KNOWN_GENES
from helpers import report, translate, START_CODONS, result_file


def check_known_genes(genome, genes, sequences):
    """Translate a few famous genes and compare their length with the table."""
    report("\nMethionine check on known genes (last 12 upstream bases + start codon):")
    rows = []
    for name in KNOWN_GENES:
        row = genes[genes["Symbol"] == name]
        if len(row) == 0:
            continue
        k = row.index[0]
        cds = genome[genes.loc[k, "Begin"] - 1: genes.loc[k, "End"]]
        protein = translate(cds)
        first_aa = "M" if cds[:3] in START_CODONS else protein[0]
        n_aa = len(cds) // 3 - 1                     # minus 1 for the stop codon
        ok = (n_aa == int(genes.loc[k, "Protein length"]))
        rows.append([name, genes.loc[k, "Locus tag"], sequences[k][88:100].lower() + " " + sequences[k][100:],
                     first_aa + protein[1:10], int(genes.loc[k, "Protein length"]), n_aa,
                     "OK" if ok else "CHECK"])
        report(f"  {name:5s} {sequences[k][88:100].lower()} {sequences[k][100:]} -> {first_aa}{protein[1:10]}...  "
               f"length {n_aa} aa (table {int(genes.loc[k, 'Protein length'])})")
    table = pd.DataFrame(rows, columns=["gene", "locus_tag", "last12_upstream_and_start_codon",
                                        "translated_N_terminus", "protein_length_table",
                                        "translated_length", "check"])
    table.to_csv(result_file("tables", "step0_methionine_check_known_genes.csv"), index=False)
    return table


def check_all_genes(genome, genes, start_codons):
    """Count the start codons and check that every CDS translates to the table length."""
    codon_counts = Counter(start_codons)
    n_met = sum(c in START_CODONS for c in start_codons)
    n_length_ok = 0
    for k in range(len(genes)):
        cds = genome[genes.loc[k, "Begin"] - 1: genes.loc[k, "End"]]
        protein = translate(cds)
        if len(cds) % 3 == 0 and len(protein) - 1 == genes.loc[k, "Protein length"] and protein.endswith("*"):
            n_length_ok += 1
        else:
            report(f"  not matching the table: {genes.loc[k, 'Locus tag']} {genes.loc[k, 'Name']}")
    report(f"Start codons of the {len(genes)} genes: {dict(codon_counts.most_common())}")
    report(f"Start codon read as Methionine: {n_met}/{len(genes)}")
    report(f"Whole CDS translates to the protein length in the table: {n_length_ok}/{len(genes)}")
    return codon_counts, n_met, n_length_ok


def run_methionine_check(d):
    """d is the dictionary from run_step0()."""
    known = check_known_genes(d["genome"], d["genes"], d["sequences"])
    codon_counts, n_met, n_length_ok = check_all_genes(d["genome"], d["genes"], d["start_codons"])
    return {"known_genes": known, "codon_counts": codon_counts, "n_met": n_met, "n_length_ok": n_length_ok}
