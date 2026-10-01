# =============================================================================
#  helpers.py - small functions that every step uses:
#  printing (and keeping a log), reading/writing FASTA files, translating DNA
#  into protein, and saving the summary of all my numbers.
# =============================================================================
import os
import json

from config import RESULTS_DIR, UPSTREAM

LOG = []          # every line I print, so I can save it as run_log.txt at the end


def report(*text):
    """Print a line and also keep it for the log file."""
    line = " ".join(str(t) for t in text)
    print(line)
    LOG.append(line)


def heading(title):
    """Print a section heading, for example ===== Q1 =====."""
    report(("\n" if LOG else "") + "=" * 70)       # blank line before every heading except the first
    report(title)
    report("=" * 70)


def clear_log():
    LOG.clear()


def make_result_folders():
    for sub in ("sequences", "tables", "figures"):
        os.makedirs(os.path.join(RESULTS_DIR, sub), exist_ok=True)


def result_file(*parts):
    """Path of a file in the results folder, e.g. result_file("tables", "x.csv")."""
    return os.path.join(RESULTS_DIR, *parts)


# ---- FASTA files ---------------------------------------------------------------
def read_fasta(path):
    """Read a FASTA file -> dictionary {sequence name: DNA string}."""
    sequences = {}
    name = None
    parts = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(parts).upper()
                name = line[1:].split()[0]      # first word of the header is the accession
                parts = []
            elif line:
                parts.append(line)
    sequences[name] = "".join(parts).upper()
    return sequences


def write_fasta(path, records, width=70):
    """records is a list of (header, sequence)."""
    with open(path, "w", encoding="utf-8") as f:
        for header, seq in records:
            f.write(">" + header + "\n")
            for k in range(0, len(seq), width):
                f.write(seq[k:k + width] + "\n")


def write_lines(path, lines):
    """Write a list of strings to a text file, one per line."""
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


# ---- DNA -> protein --------------------------------------------------------------
# standard genetic code, codons in TCAG order
CODON_TABLE = {}
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_n = 0
for b1 in "TCAG":
    for b2 in "TCAG":
        for b3 in "TCAG":
            CODON_TABLE[b1 + b2 + b3] = _aa[_n]
            _n += 1

# In bacteria these codons can start a protein. They are all read by the special
# initiator tRNA that carries (formyl-)Methionine, so the first amino acid is M.
START_CODONS = ["ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA"]


def translate(dna):
    return "".join(CODON_TABLE.get(dna[k:k + 3], "X") for k in range(0, len(dna) - 2, 3))


def upstream_position(index):
    """Index 0..99 in the 103-base sequence -> upstream position 100..1
    (1 = the base just before the start codon). I plot it as -u."""
    return UPSTREAM - index


# ---- saving at the end -------------------------------------------------------------
def save_summary(summary):
    with open(result_file("summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))


def save_log():
    with open(result_file("run_log.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")
