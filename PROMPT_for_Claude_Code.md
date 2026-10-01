Read @HANDOVER.md first. It has the full context of my BM4322 assignment, the decisions I made, the numbers my code must produce and how the code is organised. The always-on rules are in CLAUDE.md.

My code is already split into short modules in `BM4322_A1_220472T/code/`. `main.py` runs everything, and `assignment1_notebook.ipynb` runs the same functions step by step. Please help me check it on my Windows computer and understand it:

1. Check my setup. Tell me whether numpy, pandas, matplotlib, scipy and ipykernel are installed. If something is missing, give me the exact install command.
2. Run `python main.py` from `BM4322_A1_220472T/code`, then run `python compare_results.py` from the project root. It must print IDENTICAL. If it does not, show me the differences and stop.
3. Help me run the notebook in VS Code (select the Python kernel, then Run All). Afterwards run `python compare_results.py` again to confirm the notebook gives the same results.
4. Explain each code file to me in simple words: what it does, its main functions, and how the dictionaries `d`, `met`, `q1` ... `q5` carry the results from one step to the next (HANDOVER.md section 9 has the overview).
5. Give me 10 questions my lecturer might ask about my code and method in a viva, each with a short answer.

Rules:
- Do not change any method, parameter, random seed, output file names, figure names or summary keys.
- If I ask you to change the code later, keep it simple and in my voice (first-person comments, no classes, no new libraries), keep every file short, and run `compare_results.py` after every change.
- I am on Windows: use `os.path.join` for paths and `encoding="utf-8"` when writing text files.
- If something in HANDOVER.md does not match what you find in the files, stop and ask me.
