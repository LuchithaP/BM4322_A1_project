# =============================================================================
#  config.py - all file paths and settings of my assignment in one place.
#  Every other file imports its settings from here, so if I want to change a
#  number (for example the number of genes) I only change it here.
#
#  BM4322 Genomic Signal Processing - Assignment 1 - Index No. 220472T
# =============================================================================
import os

# ---- files ------------------------------------------------------------------
CODE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(CODE_DIR, "..", "data")
RESULTS_DIR = os.path.join(CODE_DIR, "..", "results")

GENOME_FILE = os.path.join(DATA_DIR, "GCA_022869965.1_ASM2286996v1_genomic.fna")
TABLE_FILE = os.path.join(DATA_DIR, "ncbi_dataset.tsv")

# The TSV uses the RefSeq name of the chromosome (NZ_CP064385.1) but the GenBank
# FASTA calls the same chromosome CP064385.1. RefSeq copies the GenBank sequence,
# so the positions in the table are valid for this FASTA (I checked that the two
# FASTA files from NCBI have the same length and the same sequence).
CHROM_IN_TABLE = "NZ_CP064385.1"
CHROM_IN_FASTA = "CP064385.1"

# ---- Step 0: which genes and how many bases ----------------------------------
N_GENES = 1000          # number of sense strand genes
UPSTREAM = 100          # bases before the start codon
DOWNSTREAM = 3          # bases from the start codon (the start codon itself)
KNOWN_GENES = ["murJ", "rpoD", "crp", "hns", "fis", "ompC", "hilA", "tolC", "trpA"]

# ---- Q1: local search ------------------------------------------------------------
QUERY = "WWWW"          # W = A or T (IUPAC code)
MATCH = 1               # local search scores - the same values as the local
MISMATCH = -1           # search code given in the module
GAP = -2

# ---- Q3 and Q4: PPM and statistical alignment --------------------------------------
WINDOW = 6              # length of the sequence used for the PPM (the Pribnow box is 6 bp)
PSEUDO = 0.01           # added to every count so that log(0) never happens
THRESHOLDS = [-1, -2, -3, -4, -5]   # statistical alignment thresholds (w.r.t. the consensus)
MAIN_T = -2             # threshold I used when comparing positions in Q5

# ---- Q5: chance-level control ----------------------------------------------------
N_SHUFFLE = 100         # shuffled copies of every sequence
SEED = 4322             # fixed random seed so the results are repeatable

# ---- letters and colours -----------------------------------------------------------
BASES = "ACGT"
WEAK = "AT"             # the bases that W stands for
BLUE, ORANGE, GREEN, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#52514e"
