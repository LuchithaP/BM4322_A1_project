# Handover: BM4322 Assignment 1 – Promoter Discovery in Bacteria

This document gives Claude Code (and me) the full context of my assignment so the work can continue
in VS Code. The code is already split into short modules, with `main.py` to run everything and a
Jupyter notebook that runs the same functions step by step. **Every result is identical to the
numbers in my report**: `summary.json`, all tables, all sequences and all 7 figures were checked
byte by byte against the earlier single-script version.

---

## 1. How to use this document

1. Unzip the project bundle and open the folder `BM4322_A1_project` in VS Code (File → Open Folder).
   Windows "Extract All" may create `BM4322_A1_project\BM4322_A1_project`; open the inner folder
   (the one that contains `CLAUDE.md`).
2. Open any file, then open Claude Code (Spark icon at the top right of the editor).
3. Paste the prompt from `PROMPT_for_Claude_Code.md`. It points Claude Code to this file.
4. `CLAUDE.md` holds the short, always-on rules. Claude Code loads it automatically because it is
   in the folder you opened.

---

## 2. Project at a glance

| Item | Value |
|---|---|
| Module | BM4322 Genomic Signal Processing (University of Moratuwa), Semester 7, 2022 batch |
| Assignment | Assignment 1, individual, 50% of the final grade, due **2026-10-02** |
| Student | Index No. **220472T** (name not in the files yet: `\studentname{}` in `report_latex/main.tex`) |
| Organism | *Salmonella enterica* subsp. *enterica* strain PartC-Senterica-RM8376 |
| Accession | **ASM2286996v1** (GenBank GCA_022869965.1, RefSeq GCF_022869965.1) |
| Chromosome | GenBank `CP064385.1` = RefSeq `NZ_CP064385.1`, 4,857,492 bp, GC 52.22% |
| Language | Python 3 (numpy, pandas, matplotlib, scipy; ipykernel for the notebook). Report in LaTeX. |
| Deliverables | 1) report PDF (LMS); 2) a separate .zip with the .fasta, .tsv, derived sequences and code |

---

## 3. The assignment brief (verbatim)

> Salmonella enterica is a rod-shaped, flagellate, facultative anaerobic, Gram-negative specie of
> bacteria commonly found in the gut of poultry. It causes Salmonellosis in humans. For the
> assignment the GenBank genome (.fasta) and its protein table (.tsv) of the given accession are
> required. Obtain 100 bases upstream and 3 bases downstream for 1000 genes of the sense strand.
> Initially for a few known genes use the presence of Methionine at the beginning of each coding
> sequence to verify your code. For the obtained sequence perform the following operations.
>
> 1. Perform a standard local search (for an intact query) to locate the WWWW promoter within each
>    sequence. Obtain the percentage of genes with potential promoters and the distribution of the
>    upstream position.
> 2. If a WWWW promoter is found find the number of consecutive Ws. Obtain the distribution of the
>    number consecutive Ws.
> 3. Obtain the sequence starting from each upstream position of (1) which is 6 bases long
>    (regardless of its C/G content). Calculate the position probability matrix.
> 4. Using the position probability matrix of (3) obtain a statistical alignment of the same
>    sequences.
> 5. Compare the percentage of promoters detected and position distribution results of the two
>    searches. Comment on the results.
>
> Deliverables: 1. Brief report (.pdf via LMS) with results and discussion of the five questions.
> 2. Soft copies of the .fasta, .tsv and any derivative sequences (.txt or .fasta) and code have to
> accompany the LMS submission in a separate .zip file.
>
> (My row in the table: Number 4 – Student ID 220472T – Accession ASM2286996v1.)

The original PDF is in `assignment_brief/`.

---

## 4. Folder layout of the bundle

```
BM4322_A1_project/
├── CLAUDE.md                     short rules, loaded automatically by Claude Code
├── HANDOVER.md                   this file
├── PROMPT_for_Claude_Code.md     the prompt to paste
├── compare_results.py            checks results/summary.json against reference_results/
├── reference_results/summary.json   the numbers used in the report
├── assignment_brief/             the assignment PDF
├── BM4322_A1_220472T/            <- the folder that is zipped for submission
│   ├── README.txt
│   ├── data/
│   │   ├── GCA_022869965.1_ASM2286996v1_genomic.fna    GenBank genome (FASTA)
│   │   └── ncbi_dataset.tsv                            protein table
│   ├── code/                     modules + main.py + assignment1_notebook.ipynb (section 9)
│   └── results/                  made by the code: sequences/ tables/ figures/ summary.json run_log.txt
└── report_latex/
    ├── main.tex                  the report (student voice)
    ├── figures/                  copies of results/figures/*.png used by main.tex
    └── 220472T_BM4322_Assignment1_Report.pdf          current compiled report (11 pages)
```

---

## 5. Data files and how they are used

| File | What it is | How it is used |
|---|---|---|
| `GCA_…_genomic.fna` | GenBank genome, FASTA. Two records: `CP064385.1` (chromosome) and `CP064386.1` (plasmid, not used) | The DNA "signal". All sequences are cut from the chromosome record. |
| `ncbi_dataset.tsv` | NCBI gene/protein table, 4,850 rows. Columns used: `Accession`, `Begin`, `End`, `Orientation`, `Gene Type`, `Protein length`, `Locus tag`, `Symbol`, `Name` | Gene positions (1-based), strand, gene type; protein length for the Methionine check. |

Facts that matter:
- The TSV says `NZ_CP064385.1` (RefSeq name). The GenBank FASTA says `CP064385.1`. The sequence is
  identical (checked with MD5), so the TSV coordinates apply directly to the GenBank FASTA.
- Positions in the TSV are **1-based**. Python is 0-based.
- One + strand gene crosses the origin of the circular chromosome (`Begin` > `End`). It is skipped.
- `Protein length` is read by pandas as float (NaN for RNA genes). `Symbol` can be empty → `fillna("")`.
- Gene counts: 4,850 in the table (4,612 protein-coding, 113 pseudogenes, 125 RNA genes); chromosome
  4,733 (+ strand 2,280, − strand 2,453); + strand protein-coding 2,154 (2,153 usable).
- The NCBI download also had GFF/GTF/GBFF, CDS and protein FASTA, RefSeq copy and JSON reports.
  They are **not used** (same annotation in other formats, or no upstream DNA).

---

## 6. Method and decisions (do NOT change any of these)

**Position convention.** In the 103-base sequence, index `i = 0..99` is upstream and `i = 100..102`
is the start codon. Upstream position `u = 100 - i` (u = 1 is the base just before the start codon).
Reports and plots show it as `-u`. The position of a hit = the position of its first (5') base.

**Step 0 – genes and sequences**
- Genes: TSV rows with `Accession == "NZ_CP064385.1"`, `Orientation == "plus"`,
  `Gene Type == "protein-coding"`, `Begin < End`; sort by `Begin` (stable); take the first 1000
  (murJ, Begin 571 … acpT, End 2,496,884).
- Sequence: `genome[Begin-101 : Begin+2]` (103 bases). `upstream = seq[:100]`, `start_codon = seq[100:]`.
- Methionine check: known genes `murJ, rpoD, crp, hns, fis, ompC, hilA, tolC, trpA` (all ATG → M,
  translated length = TSV protein length). All 1000: start codon in
  `["ATG","GTG","TTG","CTG","ATT","ATC","ATA"]` counts as Met; whole CDS `genome[Begin-1:End]`
  translated with the standard code; OK if `len(protein)-1 == Protein length` and ends with `*`.

**Q1 – standard local search (Smith–Waterman)**
- Query `WWWW` (W = A or T). Score +1 if base in "AT", else −1. Linear gap −2 (module values).
- Matrix `H` (101 × 5), `H[i,j] = max(0, diag + s, up + GAP, left + GAP)`.
- `best = H.max()`; start of the alignment by traceback from `np.argmax(H)` = **first** maximum
  → the 5'-most hit. Intact if `best == 4`.
- Search only the 100 upstream bases. Cross-check: intact hit == first match of regex `[AT]{4}`
  (an `assert` in the code).

**Q2 – consecutive Ws:** count A/T bases from the hit start (within the 100 upstream bases).
Compare with the geometric distribution `P(L=n) = p^(n-4) (1-p)`, `p` = A/T fraction of all
upstream bases (0.54913).

**Q3 – PPM:** hexamer = `sequences[k][start:start+6]` for each intact hit. Counts 4×6 (rows A,C,G,T).
`PPM = (counts + 0.01) / (N + 0.04)`, `log_ppm = ln(PPM)`, consensus = argmax per column,
`consensus_score = sum of column max of log_ppm`, information `I_j = Σ p log2(p/0.25)`.

**Q4 – statistical alignment:** 95 windows over the 100 upstream bases (−100 … −6).
`S = Σ_j log_ppm[base_j, j]`; best window = `np.argmax` (first max on ties);
`rel = S_best − consensus_score`; detected if `rel >= T` for `T in [-1,-2,-3,-4,-5]`; `MAIN_T = -2`.

**Q5 – comparison and controls**
- Same position / within 2 bases / shift; `scipy.stats.ks_2samp(q1_positions, q4_positions)`.
- Shuffled control: `random.seed(4322)` **immediately before** the loop
  `for seq in upstream: for _ in range(100): "".join(random.sample(seq, len(seq)))`
  (in `q5_comparison.shuffle_sequences`). The order of these random calls must not change.
  WWWW in shuffled sequences is found with the regex (SW is too slow; equivalence is asserted on the
  real data).
- Hits by region (−100…−26, −25…−14, −13…−9, −8…−6) for real vs shuffled.
- Overlap of each upstream window with any annotated chromosome gene (coverage mask, origin-crossing
  gene included); same-strand neighbour within 100 bases = "operon-like".
- Detection vs threshold sweep 0 to −12 in steps of −0.25 (real and shuffled).
- Chance expectation: `97 × p^4` WWWW per sequence; P(at least one) by a small Markov-chain loop.

---

## 7. Outputs (names must stay the same – the report uses them)

`results/sequences/`: `upstream100_downstream3_1000genes.fasta`, `upstream100_downstream3_1000genes.txt`,
`Q3_hexamers.fasta`, `Q3_hexamers.txt`, `Q4_best_6mers.fasta`

`results/tables/`: `gene_list_1000.csv`, `step0_methionine_check_known_genes.csv`,
`Q1_local_search_results.csv`, `Q1_position_distribution.csv`, `Q1_genes_without_WWWW.csv`,
`Q2_consecutive_W_distribution.csv`, `Q3_frequency_matrix.csv`, `Q3_PPM.csv`,
`Q4_statistical_alignment_results.csv`, `Q4_Q5_detection_vs_threshold.csv`, `Q5_hits_by_region.csv`

`results/figures/` (used by `main.tex`): `fig_Q1_positions.png`, `fig_Q1_all_occurrences.png`,
`fig_Q2_consecutive_W.png`, `fig_Q3_PPM_logo.png`, `fig_Q4_positions.png`, `fig_Q5_threshold.png`,
`fig_Q5_positions.png`

`results/summary.json` (81 keys) and `results/run_log.txt`.

---

## 8. Key results (the code must give exactly these)

| What | Value |
|---|---|
| Start codons (1000 genes) | ATG 912, GTG 66, TTG 18, ATT 2, CTG 2; Met 1000/1000; CDS length OK 999/1000 (IUJ38_RS08235, IS3 transposase, frameshift) |
| Upstream A/T | 54.913% (genome 47.78%) |
| Q1 intact WWWW | 995/1000 = 99.5% (best score 4: 995, 3: 5) |
| Q1 position | median −90, IQR −97 to −77, 135 hits at −100, 720 at −80 or further |
| WWWW per sequence | mean 12.224, median 11, range 0–42, total 12,224; last-hit median −14 |
| Genes without WWWW | holB, IUJ38_RS03200, pduS, srlE, gltB |
| Q2 runs | mean 5.489, median 5, range 4–22; {4:381, 5:247, 6:157, 7:88, 8:53, 9:33, 10:16, 11:6, 12:4, 13:6, 14:2, 22:2}; ≥8: 122 vs 90.47 expected |
| Q3 | 995 hexamers, 239 distinct; consensus ATTTTT; score −5.0122; information 1.001, 1.007, 1.000, 1.000, 0.047, 0.009 (total 4.064 bits) |
| Q3 counts | A [517,449,483,485,293,268]; C [0,0,0,0,220,251]; G [0,0,0,0,161,203]; T [478,546,512,510,321,273] |
| Q4 detection | T=−1: 990 (99.0%); T=−2…−5: 995 (99.5%) |
| Q4 position | median −51, IQR −76 to −28; top windows ATTTTT 113, ATTATT 64, ATTTTA 56, TTTTTT 43; all-A/T 682 |
| Q4 other | no-hit scores −10.90 … −11.19; worst WWWW window −1.369; cost of C at position 1 −10.853 |
| Q5 | same position 179, within 2 bases 241, median shift 30, KS D = 0.5317 (p ≈ 5.1e−129); run at stat hit 6.384 vs 5.489 |
| Shuffled control | WWWW 97.784%; 9.698 per sequence; statistical detection T=−1 96.233%, T≤−2 97.638%; medians −88 (Q1) and −54 (Q4) |
| Regions (stat real vs shuffled) | −100…−26 77.59/79.11; −25…−14 16.08/12.53; −13…−9 1.61/5.13; −8…−6 4.72/3.23 |
| Overlap with other genes | 466 genes (428 same-strand neighbour within 100 bases); WWWW in 99.14% of them vs 99.81% of free ones |
| Chance expectation | 8.82 WWWW per 100 bases; P(at least one) 0.9936 |

---

## 9. Code structure (done)

All files are in `BM4322_A1_220472T/code/`. Each step is a `run_...()` function that takes the
results of the earlier steps (dictionaries) and returns its own dictionary. `main.py` and the
notebook call exactly the same functions in the same order, so they give the same results.

| File | Lines (about) | Main functions | Returns |
|---|---|---|---|
| `config.py` | 50 | constants only (paths, N_GENES, UPSTREAM, QUERY, scores, WINDOW, PSEUDO, THRESHOLDS, MAIN_T, N_SHUFFLE, SEED, colours) | – |
| `helpers.py` | 111 | `report`, `heading`, `clear_log`, `make_result_folders`, `result_file`, `read_fasta`, `write_fasta`, `write_lines`, `translate`, `upstream_position`, `save_summary`, `save_log` | – |
| `step0_data.py` | 108 | `load_genome`, `load_table`, `choose_genes`, `cut_sequences`, `gene_coverage`, `save_sequences`, **`run_step0()`** | `d`: genome, gc, table, on_chrom, genes, sequences, upstream, start_codons, coverage, p_w |
| `step0_methionine_check.py` | 62 | `check_known_genes`, `check_all_genes`, **`run_methionine_check(d)`** | `met`: known_genes, codon_counts, n_met, n_length_ok |
| `q1_local_search.py` | 123 | `score`, `local_search`, `search_all`, `genes_without_hit`, **`run_q1(d)`** | `q1`: sw_score, hit_start, intact, found, q1_pos, q1_counts, all_hits, n_occ, occ_counts, last_pos, no_hit_rows |
| `q2_consecutive_w.py` | 53 | `count_consecutive_w`, **`run_q2(d, q1)`** | `q2`: runs, run_counts, expected, obs_ge8, exp_ge8, obs_ge15, exp_ge15, longest |
| `q3_ppm.py` | 71 | `get_hexamers`, `make_ppm`, `information_content`, **`run_q3(d, q1)`** | `q3`: hexamers, counts, ppm, log_ppm, consensus, consensus_score, info |
| `q4_statistical_alignment.py` | 85 | `to_numbers`, `stat_align(seq, log_ppm)`, `worst_wwww_score`, **`run_q4(d, q1, q3)`** | `q4`: best_index, best_rel, stat_pos, best_hex, detected, det_pos, top_best, all_w, worst_wwww, no_hit_scores |
| `q5_comparison.py` | 160 | `compare_positions`, `shuffle_sequences`, `shuffled_control`, `hits_by_region`, `overlap_with_other_genes`, `threshold_sweep`, `chance_expectation`, **`run_q5(d, q1, q2, q3, q4)`** | `q5`: both, same, close, shift, ks, run_at_stat, control, region_rows, inside, same_strand_close, found_inside, found_free, sweep, real_curve, shuf_curve, expected_hits, p_at_least_one |
| `plots.py` | ~230 | `q1_positions`, `all_occurrences`, `q2_runs`, `draw_logo`, `q3_ppm_logo`, `q4_positions`, `q5_threshold`, `q5_positions`, **`save_all_figures(...)`** | figures (each function returns a figure) |
| `summary.py` | 71 | **`build_summary(d, met, q1, q2, q3, q4, q5)`** | dict with the 81 keys of `summary.json` |
| `main.py` | 46 | runs Step 0 → Q5, saves figures, summary and log | – |
| `assignment1_notebook.ipynb` | 29 cells | same calls as `main.py`, with explanations, tables, a demo of `local_search` and `stat_align` on one sequence, and the figures under each step | – |

Data flow: `run_step0` → `d`; `run_methionine_check(d)` → `met`; `run_q1(d)` → `q1`;
`run_q2(d, q1)`; `run_q3(d, q1)`; `run_q4(d, q1, q3)`; `run_q5(d, q1, q2, q3, q4)`;
`plots.save_all_figures(d, q1, q2, q3, q4, q5)`; `build_summary(d, met, q1, …, q5)` → `summary.json`.

---

## 10. How to check that nothing changed

1. Run `python main.py` from `BM4322_A1_220472T/code` (or Run All in the notebook).
2. From the project root run `python compare_results.py`. It compares
   `BM4322_A1_220472T/results/summary.json` with `reference_results/summary.json` and must print
   `IDENTICAL` (tiny float noise < 1e-9 is ignored).
3. Check that `results/` has exactly the file names in section 7 and that the figures look the same.

---

## 11. If the code is changed later

- Run section 10 after every change.
- `BM4322_A1_220472T/README.txt`: keep the file list and run commands up to date (student voice).
- `report_latex/main.tex`: Appendix A lists the code files; Appendix B shows `score()` +
  `local_search()` from `q1_local_search.py` and `stat_align()` + the best-window loop from
  `q4_statistical_alignment.py` – update them if those functions change. Copy
  `results/figures/*.png` into `report_latex/figures/`. Compile with `latexmk -pdf main.tex` (or on
  Overleaf).
- Rebuild the submission zip from the project root after deleting every `__pycache__` folder:

```
python -c "import shutil; shutil.make_archive('220472T_BM4322_Assignment1_Code_and_Data', 'zip', root_dir='.', base_dir='BM4322_A1_220472T')"
```

---

## 12. Environment (my computer)

- Windows, VS Code, Claude Code extension. Terminal is PowerShell.
- Python 3.10 or newer. Install: `python -m pip install numpy pandas matplotlib scipy ipykernel`
  (use `py` instead of `python` if `python` is not found).
- `main.py` takes about 10 seconds; the notebook about 15 seconds.
- In VS Code, open the notebook from the `code` folder and pick the Python interpreter as kernel.

---

## 13. Style guide: "student perspective"

The report and the code must read as if I wrote them as a final-year student.
- First person and plain English: "I read the table…", "I skip this gene because…".
- Comments explain **why**, briefly. No marketing words (robust, leverage, comprehensive, seamless,
  delve, crucial). No emojis.
- Simple code: functions and loops, no classes, no argparse, no logging module, no type-hint
  clutter, no new libraries (only numpy, pandas, matplotlib, scipy and the standard library).
- Short docstrings (1–2 lines). Each module starts with a short comment saying which question it
  answers.
- The report (LaTeX) uses "I", simple sentences, numbers from `summary.json`, positions written
  as `$-90$`, sequences with `\seqs{...}`, gene names with `\gene{...}`.

---

## 14. Known pitfalls

- 1-based (TSV) vs 0-based (Python): the slice is `genome[Begin-101 : Begin+2]`.
- `NZ_CP064385.1` (TSV) vs `CP064385.1` (GenBank FASTA) – same chromosome.
- The origin-crossing gene (`Begin > End`) is excluded from the 1000 genes but **included** in the
  coverage mask used for the overlap analysis.
- `np.argmax` returns the first maximum: this defines the "first hit" rule in Q1 and the tie rule
  in Q4. Do not replace it with something that picks a different tie.
- `random.seed(SEED)` is called inside `shuffle_sequences()` right before the loop; nothing else may
  use `random` there.
- Only `main.py` calls `matplotlib.use("Agg")`. If a module did it, the notebook could not show figures.
- In the notebook, after editing a `.py` file, restart the kernel (or the old version stays loaded).
- `json.dump` uses numpy values – keep the converter in `helpers.save_summary`.
- Write text files with `encoding="utf-8"` (Windows default is not UTF-8).

---

## 15. Glossary (for explaining the code to me)

- **Promoter**: DNA in front of a gene where RNA polymerase binds. In bacteria: −35 box `TTGACA`
  and −10 (Pribnow) box `TATAAT`, about 17 bases apart, about 10 bases before transcription starts.
- **W**: IUPAC code for A or T ("weak", 2 hydrogen bonds). `TATAAT` = `WWWWWW`.
- **Start codon / Methionine**: genes start with ATG (sometimes GTG/TTG); the first amino acid is Met.
- **Shine–Dalgarno (SD)**: `AGGAGG`, ~8–13 bases before the start codon, where the ribosome binds.
- **Smith–Waterman**: local alignment by dynamic programming; the 0 in the recurrence lets the match
  start anywhere. "Intact query" = whole query matched with no gap or mismatch (score 4 here).
- **PFM / PPM**: counts / probabilities of each base at each position of the aligned 6-mers.
- **Pseudocount**: small number added to counts so `ln(0)` never happens.
- **Information content**: how far a position is from random, in bits (0 = random, 2 = fixed base).
- **Statistical alignment**: slide the PPM along the sequence, add log probabilities (like a matched
  filter), compare the best score with the consensus score using a threshold.
- **Shuffled control**: same bases in random order; shows what the search finds by chance.
- **KS test**: tests whether two distributions (here, hit positions) are different.

---

## 16. History

1. v1: single script + Word/PDF report (more formal voice).
2. v2: code and report rewritten in my voice; report moved to LaTeX; only the GenBank FASTA and the
   TSV are used. All report numbers were cross-checked against `summary.json`, and an independent
   re-implementation gave the same core numbers.
3. v3 (current): code split into modules + `main.py` + `assignment1_notebook.ipynb`. `summary.json`,
   every table, every sequence file and all 7 figures are byte-identical to v2. README and the
   report appendices were updated (report is 11 pages).

---

## 17. Next tasks for Claude Code

1. Check the setup on my Windows computer and run `main.py` + `compare_results.py` (must be
   `IDENTICAL`).
2. Run the notebook in VS Code and check again.
3. Explain each file to me in simple words and give me likely viva questions with short answers.
4. Only if I ask: small changes, always followed by section 10 and section 11.
