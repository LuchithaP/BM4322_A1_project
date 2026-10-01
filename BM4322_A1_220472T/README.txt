BM4322 Genomic Signal Processing - Assignment 1: Promoter Discovery in Bacteria
Index No. 220472T   |   Accession ASM2286996v1
Salmonella enterica subsp. enterica strain PartC-Senterica-RM8376
==============================================================================

This folder has my data, code and results for the assignment.

data/
  GCA_022869965.1_ASM2286996v1_genomic.fna   the GenBank genome (.fasta) that I downloaded from NCBI
  ncbi_dataset.tsv                           the protein table (.tsv) for the same genome

code/  (one file for each step of the assignment)
  main.py                        runs the whole assignment in order (Step 0 -> Q5) from the command line
  assignment1_notebook.ipynb     runs the same functions step by step, with my explanations,
                                 the results and the figures under each step
  config.py                      file paths and all the settings (numbers) I used
  helpers.py                     printing + log, reading/writing FASTA files, translating DNA
  step0_data.py                  Step 0: choose the 1000 genes and cut the 103-base sequences
  step0_methionine_check.py      Step 0: check that the start codon gives Methionine
  q1_local_search.py             Q1: Smith-Waterman local search for the intact query WWWW
  q2_consecutive_w.py            Q2: number of consecutive Ws at each hit
  q3_ppm.py                      Q3: position probability matrix (PPM)
  q4_statistical_alignment.py    Q4: statistical alignment with the PPM
  q5_comparison.py               Q5: comparison of the two searches + shuffled control
  plots.py                       the figures of my report
  summary.py                     collects all my numbers into results/summary.json

results/  (made by main.py or by the notebook - both give exactly the same files)
  sequences/upstream100_downstream3_1000genes.fasta   the 1000 sequences I cut out
                                                     (100 bases upstream + 3 bases of the start codon)
  sequences/upstream100_downstream3_1000genes.txt     same sequences, one per line
                                                     (upstream bases in capitals, start codon in small letters)
  sequences/Q3_hexamers.fasta / .txt                  the 6 bases starting at every WWWW hit (used for the PPM)
  sequences/Q4_best_6mers.fasta                       the best window of every sequence in the statistical alignment
  tables/*.csv                                        the results of every question
  figures/*.png                                       all the figures in my report
  summary.json                                        all the numbers I used in the report
  run_log.txt                                         everything the code prints

How to run
----------
I used Python 3 with numpy, pandas, matplotlib and scipy (and ipykernel to run the notebook):
    pip install numpy pandas matplotlib scipy ipykernel

Option 1 - the whole assignment in one go:
    cd code
    python main.py
Option 2 - step by step:
    open code/assignment1_notebook.ipynb in VS Code or Jupyter and choose "Run All"
    (the notebook has to be opened from the code folder so it can import my .py files)

Both take about 10-15 seconds and write everything into the results folder.

Notes
-----
- Positions are given relative to the start codon: -1 is the base just before the start codon
  and -100 is the farthest base. The position of a hit is the position of its first base.
- The protein table calls the chromosome NZ_CP064385.1 (RefSeq name) and the GenBank FASTA calls it
  CP064385.1. They are the same sequence, so the table positions work with the GenBank genome.
- The shuffled control uses a fixed random seed (4322), so the results are the same every time.
