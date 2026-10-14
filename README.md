# Certifying zero false alarms: code and reproduction

Code accompanying *When can a finite screen certify zero false alarms? Distribution-free
impossibility, the detectability gap, and a validation budget for off-target screening*, by
Amir Daneshpajouh, Megan Fowler and Kay C. Wiese (School of Computing Science, Simon Fraser
University).

The paper is a theory paper. This repository computes every number the paper reports: the
exact sample sizes, the validation-budget table and figure, and the conceptual figure, and it
checks the sufficiency construction of Theorem 2(i) by exact computation. The
one input that is not computed here from first principles is the enumeration count N for the
worked example, which comes from scanning a reference genome; `enumerate_sites.py` does that
scan, `data/exacel_enumeration_summary.json` records its output, and `docs/ENUMERATION.md`
documents it.

## What is here

| Path | What it does |
|---|---|
| `m0_search.py` | Exact integer sample-size search. The authority for every sample size in the paper. |
| `test_m0_search.py` | 12 tests, including the exact boundary cases. |
| `boundary_exact_test.py` | Exact-rational check of the sample-size formula against the Neyman-Pearson umbrella criterion, at the boundary cases where floating point fails. |
| `occupancy_curve_test.py` | Compares two closed forms with the exact finite-population requirement. |
| `test_certificate.py` | Exact rational check of Theorem 2(i) on random and adversarial discrete classes: the tolerance statement, the type-I bound, and the power bound (3 tests, about a minute). |
| `enumerate_sites.py` | Scans a reference genome for candidate off-target windows of a spacer. |
| `test_enumerate_sites.py` | Tests the scan against an independent naive implementation on planted sequences. |
| `data/exacel_enumeration_summary.json` | Output of the scan for the exa-cel spacer on GRCh38: window counts by mismatch, distinct 23-mers, the on-target window, and the size and SHA-256 of every reference contig file. |
| `data/exacel_windows.tsv.gz` | The 28,176 windows: contig, strand, 1-based forward start, 23-mer, mismatches. |
| `gamma_from_enumeration.py` | Derives the detectability-gap bound k/N and the budgets from that summary. |
| `make_budget_curve.py` | Produces Figure 2 and Table 1 of the paper (`budget_curve.pdf`, `table_main.tex`, `table_appendix.tex`, `budget_table.json`). |
| `make_schematic.py` | Produces Figure 1 of the paper (`contamination.pdf`). |
| `docs/ENUMERATION.md` | How N is obtained, the unit conventions, and the counts. |
| `docs/GAMMA_EXACEL.txt` | Printed output of `gamma_from_enumeration.py`. |

Python 3.9 or later. `numpy` is needed for the scan and `matplotlib` for the figures;
everything else uses the standard library.

## Reproducing the paper's numbers

```
python3 test_m0_search.py          # 12 tests, all must pass
python3 test_enumerate_sites.py    # 3 tests, all must pass
python3 test_certificate.py        # 3 tests, about a minute
python3 make_budget_curve.py       # Figure 2, Table 1 and the appendix table
python3 make_schematic.py          # Figure 1
python3 gamma_from_enumeration.py  # the printed tables behind Table 1
```

To redo the enumeration from a reference (about four minutes and several GB of memory for
the whole of GRCh38, one file per contig in a directory):

```
python3 enumerate_sites.py --ref-dir DIR_WITH_CONTIG_FILES --out-json summary.json
```

## The one thing worth knowing before reading the code

The least admissible sample size is `m0 = min{ m : (1 - gamma)^m <= delta1 }`. As a
statement about real numbers that equals `ceil(log(delta1) / log(1 - gamma))`, and that is
how it is usually printed. **Evaluating it in floating point is wrong at exact boundaries.**
When `delta1 = (1 - gamma)^k` the ceiling can come out one too large; the smallest case is
`gamma = 1/5, delta1 = 0.64`, where the answer is 2 and the float formula returns 3 when delta1 is entered as the double 0.64 (entered as the double `(4/5)**2` it returns 2; the result depends on how the inputs are represented, which is why the code takes rationals). A scan
of 36 exact boundary cases found 9 that overshoot.

`m0_search.py` therefore takes rationals and decides the comparison exactly, using
arbitrary-precision decimal arithmetic with a verified error bound and falling through to
exact integer comparison on any tie or near-tie. Floating point is used only to guess a
starting bracket, which is then verified.

The same applies to the finite-population budget. The closed form `1 - delta1**(1/k)` is
the leading-order expression for the hypergeometric zero-defect requirement with the
continuity correction dropped, and it overshoots the exact threshold by up to two sites on
the candidate set used in the paper. The paper reports the exact integer from
`m0_hypergeom`.

## Scope

The budgets are validation requirements for a *statistical* certificate. They are not
evidence about the clinical safety of any therapy. The candidate-set enumeration is
reference-based, so a site absent from the reference, including a site created by a common
variant, is outside the frame the argument bounds, and the assay is assumed to detect
cleavage at any validated site where it occurs.

The sufficiency procedure implemented here is classical. The sample-size requirement is the
one-sided distribution-free tolerance bound of Wilks (1941), the same quantity appears as
the minimum class-0 sample size in the Neyman-Pearson umbrella algorithm of Tong, Feng and
Li (2018), and the finite-population form belongs to the surveillance literature (Cameron
and Baldock 1998; Johnson et al. 2004). The paper's contribution is what the procedure
certifies and on which classes, not the procedure.

## Licence

Code: MIT licence (see `LICENSE`). Data, tables and figures: CC BY 4.0 (see `data/LICENSE-DATA.md`).
