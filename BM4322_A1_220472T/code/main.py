# =============================================================================
#  BM4322 Genomic Signal Processing - Assignment 1
#  Promoter Discovery in Bacteria
#
#  Index No.  : 220472T
#  Accession  : ASM2286996v1  (GenBank GCA_022869965.1)
#  Organism   : Salmonella enterica subsp. enterica strain PartC-Senterica-RM8376
#
#  main.py runs the whole assignment from the command line, in the same order as
#  the questions. The same functions are run step by step in the notebook
#  assignment1_notebook.ipynb.
#
#  How to run (from this folder):
#    pip install numpy pandas matplotlib scipy
#    python main.py
#  All results go to the ../results folder (sequences, tables, figures, summary).
# =============================================================================
import matplotlib
matplotlib.use("Agg")                      # only save figures, do not open windows

import helpers
import plots
from step0_data import run_step0
from step0_methionine_check import run_methionine_check
from q1_local_search import run_q1
from q2_consecutive_w import run_q2
from q3_ppm import run_q3
from q4_statistical_alignment import run_q4
from q5_comparison import run_q5
from summary import build_summary

helpers.clear_log()
helpers.make_result_folders()

d = run_step0()                    # Step 0: genes and 103-base sequences
met = run_methionine_check(d)      # Step 0: Methionine check
q1 = run_q1(d)                     # Q1: local search for WWWW
q2 = run_q2(d, q1)                 # Q2: consecutive Ws
q3 = run_q3(d, q1)                 # Q3: position probability matrix
q4 = run_q4(d, q1, q3)             # Q4: statistical alignment
q5 = run_q5(d, q1, q2, q3, q4)     # Q5: comparison + shuffled control

plots.save_all_figures(d, q1, q2, q3, q4, q5)
helpers.save_summary(build_summary(d, met, q1, q2, q3, q4, q5))
helpers.save_log()
helpers.report("\nDone - results saved in the 'results' folder.")
