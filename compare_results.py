# compare_results.py - checks that my code still gives exactly the same numbers as in my report.
# Run it from the project folder after running main.py (or the notebook):
#     python compare_results.py
# It compares BM4322_A1_220472T/results/summary.json with reference_results/summary.json
# (the numbers used in the report). Tiny rounding noise (< 1e-9) is ignored.
import json
import math


def same(a, b, tol=1e-9):
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(same(a[k], b[k], tol) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(same(x, y, tol) for x, y in zip(a, b))
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return math.isclose(a, b, rel_tol=tol, abs_tol=tol)
    return a == b


reference = json.load(open("reference_results/summary.json", encoding="utf-8"))
now = json.load(open("BM4322_A1_220472T/results/summary.json", encoding="utf-8"))
bad = [k for k in sorted(set(reference) | set(now)) if not same(reference.get(k), now.get(k))]
print("IDENTICAL" if not bad else "DIFFERENT: " + ", ".join(bad))
