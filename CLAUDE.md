# BM4322 Assignment 1 (Promoter Discovery in Bacteria) – rules for Claude Code

Student project, index 220472T, accession ASM2286996v1 (Salmonella enterica). Full context, the
expected numbers and the code structure are in `HANDOVER.md` – read it before editing anything.

## Layout
- `BM4322_A1_220472T/` is the submission folder: `data/` (GenBank .fna + .tsv), `code/`
  (modules, `main.py`, `assignment1_notebook.ipynb`), `results/` (made by the code), `README.txt`.
- `reference_results/summary.json` = the numbers in the report. `compare_results.py` checks them.
- `report_latex/main.tex` is the report; it uses the figures in `report_latex/figures/`.

## Always
- Do not change any method, parameter, random seed (4322), the order of random calls, output file
  names, figure names, or the keys of `results/summary.json`.
- After any code change: run `python main.py` in `BM4322_A1_220472T/code`, then
  `python compare_results.py` in the project root. It must print `IDENTICAL`.
- Write code and comments in the student's voice: first person, plain English, simple functions,
  no classes, no argparse, no new libraries (numpy, pandas, matplotlib, scipy only).
- Keep every code file short. Modules must not call `matplotlib.use(...)` (it breaks the notebook);
  only `main.py` does.
- Windows: build paths with `os.path.join`; write text files with `encoding="utf-8"`.
- Never edit the files in `data/`. If HANDOVER.md and the files disagree, stop and ask.

## Commands
- Install: `python -m pip install numpy pandas matplotlib scipy ipykernel`
- Run everything: `cd BM4322_A1_220472T/code` then `python main.py`
- Check numbers: `python compare_results.py` (from the project root)
- Notebook: open `BM4322_A1_220472T/code/assignment1_notebook.ipynb`, select the Python kernel, Run All
- Report: in `report_latex/`, `latexmk -pdf main.tex` (or compile on Overleaf)
