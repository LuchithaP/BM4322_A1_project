# =============================================================================
#  step0_data.py - STEP 0 of the assignment.
#  I read the genome (.fna) and the protein table (.tsv), choose the first 1000
#  protein-coding genes on the sense (+) strand and cut 100 bases upstream +
#  3 bases downstream (the start codon) of each gene.
# =============================================================================
import numpy as np
import pandas as pd

from config import (GENOME_FILE, TABLE_FILE, CHROM_IN_TABLE, CHROM_IN_FASTA,
                    N_GENES, UPSTREAM, DOWNSTREAM)
from helpers import report, heading, read_fasta, write_fasta, write_lines, result_file


def load_genome():
    """Read the GenBank FASTA and return the chromosome and its GC fraction."""
    genome = read_fasta(GENOME_FILE)[CHROM_IN_FASTA]
    gc = (genome.count("G") + genome.count("C")) / len(genome)
    report(f"Genome {CHROM_IN_FASTA}: {len(genome):,} bp, GC content {100 * gc:.2f}%")
    return genome, gc


def load_table():
    """Read the protein table and keep the genes that are on the chromosome."""
    table = pd.read_csv(TABLE_FILE, sep="\t")
    table["Symbol"] = table["Symbol"].fillna("")          # some genes have no short name
    report(f"Protein table: {len(table)} genes, types: "
           f"{ {k: int(v) for k, v in table['Gene Type'].value_counts().items()} }")
    on_chrom = table[table["Accession"] == CHROM_IN_TABLE]
    report(f"On the chromosome: {len(on_chrom)} genes "
           f"(+ strand {sum(on_chrom['Orientation'] == 'plus')}, "
           f"- strand {sum(on_chrom['Orientation'] == 'minus')})")
    return table, on_chrom


def choose_genes(on_chrom):
    """First 1000 protein-coding genes on the + strand, in the order of the genome.
    One gene crosses the end of the circular chromosome (Begin > End), I skip it."""
    genes = on_chrom[(on_chrom["Orientation"] == "plus")
                     & (on_chrom["Gene Type"] == "protein-coding")
                     & (on_chrom["Begin"] < on_chrom["End"])]
    report(f"Protein-coding genes on the + strand: {len(genes)}")
    genes = genes.sort_values("Begin", kind="stable").head(N_GENES).reset_index(drop=True)
    report(f"I use the first {N_GENES}: {genes['Symbol'].iloc[0]} (starts at {genes['Begin'].iloc[0]}) "
           f"to {genes['Symbol'].iloc[-1]} (ends at {genes['End'].iloc[-1]})")
    return genes


def cut_sequences(genome, genes):
    """The table gives 1-based positions (the first base of the genome is 1) but Python
    strings start at 0. The start codon is at Begin..Begin+2, so in Python it is
    genome[Begin-1 : Begin+2], and the 100 bases before it start at Begin-101."""
    sequences = []
    for begin in genes["Begin"]:
        seq = genome[begin - 1 - UPSTREAM: begin - 1 + DOWNSTREAM]
        assert len(seq) == UPSTREAM + DOWNSTREAM
        sequences.append(seq)
    upstream = [s[:UPSTREAM] for s in sequences]      # the 100 bases where I search for promoters
    start_codons = [s[UPSTREAM:] for s in sequences]   # the last 3 bases
    return sequences, upstream, start_codons


def gene_coverage(on_chrom, genome_length):
    """True for every genome position that is inside an annotated gene (any strand).
    I use it to see whether an upstream region is free DNA or inside another gene."""
    cover = np.zeros(genome_length + 2, dtype=bool)
    for b, e in zip(on_chrom["Begin"], on_chrom["End"]):
        if b <= e:
            cover[b:e + 1] = True
        else:                                   # the gene that crosses the origin
            cover[b:] = True
            cover[:e + 1] = True
    return cover


def save_sequences(genes, sequences, upstream, start_codons):
    records = []
    for k in range(len(genes)):
        b = genes.loc[k, "Begin"]
        header = (f"{genes.loc[k, 'Locus tag']}|{genes.loc[k, 'Symbol'] or '-'}"
                  f"|{CHROM_IN_FASTA}:{b - UPSTREAM}-{b + DOWNSTREAM - 1}(+)|start={start_codons[k]}")
        records.append((header, sequences[k]))
    write_fasta(result_file("sequences", "upstream100_downstream3_1000genes.fasta"), records, width=103)
    write_lines(result_file("sequences", "upstream100_downstream3_1000genes.txt"),
                [f"{genes.loc[k, 'Locus tag']}\t{upstream[k]}{start_codons[k].lower()}" for k in range(len(genes))])
    genes_out = genes[["Locus tag", "Symbol", "Name", "Begin", "End", "Protein length"]].copy()
    genes_out["start_codon"] = start_codons
    genes_out["sequence_103nt"] = sequences
    genes_out.to_csv(result_file("tables", "gene_list_1000.csv"), index=False)


def run_step0():
    """Do all of Step 0 and return everything the next steps need in one dictionary."""
    heading("STEP 0 - data and sequence extraction")
    genome, gc = load_genome()
    table, on_chrom = load_table()
    genes = choose_genes(on_chrom)
    sequences, upstream, start_codons = cut_sequences(genome, genes)
    save_sequences(genes, sequences, upstream, start_codons)

    # A/T content of the upstream regions (I need it later to say what is expected by chance)
    all_up = "".join(upstream)
    p_w = (all_up.count("A") + all_up.count("T")) / len(all_up)
    report(f"A/T content of the upstream regions: {100 * p_w:.1f}% (whole genome {100 * (1 - gc):.1f}%)")

    return {"genome": genome, "gc": gc, "table": table, "on_chrom": on_chrom, "genes": genes,
            "sequences": sequences, "upstream": upstream, "start_codons": start_codons,
            "coverage": gene_coverage(on_chrom, len(genome)), "p_w": p_w}
