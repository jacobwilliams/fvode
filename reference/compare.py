"""Compare FVODE step traces against the reference C CVODE traces.

For every reference trace  reference/out/trace_<name>.csv  (written by
run_reference.sh), look for the FVODE trace  reference/out/fvode_<name>.csv
(written by the matching FVODE trace program) and compare them row by row.

Exact mode (--exact; pixi task compare-exact). FVODE is built with
-ffp-contract=off, like the C reference, so a faithful translation gives
bit-for-bit identical traces. This is the fidelity check: the first differing
row points to the step where the translation departs from CVODE.

Tolerance mode (default; pixi task compare). FVODE is built with normal
optimizations (e.g. fused multiply-add), so real values drift by rounding.
* If all integer columns (steps, orders, counters, flags) match, the real
  columns must agree to a mixed tolerance (--rtol, default 1e-3):
  |a - b| <= rtol * max(|a|, |b|)  or  |a - b| <= rtol * max|column|.
* On sensitive problems a rounding difference can flip a step decision, after
  which the step sequences differ. The trace is then reported as DIVERGED
  (not a failure) if the final solution values agree to --ytol (default 1e-4,
  scaled by the column maximum) and the final step counts differ by at most
  --nsttol (default 10%).

Usage: python reference/compare.py [--exact] [--rtol R] [--ytol Y] [--nsttol N] [--strict]
  --strict : also fail when an FVODE trace is missing
"""

import argparse
import csv
import math
import sys
from pathlib import Path

OUT = Path(__file__).resolve().parent / "out"
INT_COLS = {"nst", "qu", "qcur", "nfe", "nni", "nncf", "netf", "nsetups",
            "njev", "nfels", "flag"}


def read(path: Path) -> tuple[list[str], list[list[str]]]:
    with path.open(newline="") as fh:
        rows = list(csv.reader(fh))
    return rows[0], rows[1:]


def compare(name: str, ref: Path, new: Path, args) -> bool:
    rhead, rrows = read(ref)
    nhead, nrows = read(new)
    if rhead != nhead:
        print(f"FAIL {name}: column headers differ\n  ref: {rhead}\n  new: {nhead}")
        return False

    real_cols = [k for k, c in enumerate(rhead) if c not in INT_COLS]
    colmax = {k: max((abs(float(r[k])) for r in rrows), default=0.0) for k in real_cols}

    max_diff = 0.0
    for i, (r, n) in enumerate(zip(rrows, nrows), start=1):
        for k, (col, rv, nv) in enumerate(zip(rhead, r, n)):
            if col in INT_COLS:
                if int(rv) != int(nv):
                    if args.exact:
                        print(f"FAIL {name}: row {i} (nst={r[0]}) column '{col}': "
                              f"ref={rv} fvode={nv}")
                        return False
                    return diverged(name, rhead, rrows, nrows, real_cols, colmax, i, col, args)
                continue
            a, b = float(rv), float(nv)
            if a == b:
                continue
            if args.exact:
                print(f"FAIL {name}: row {i} (nst={r[0]}) column '{col}': ref={rv} fvode={nv}")
                return False
            rel = abs(a - b) / max(abs(a), abs(b))
            scaled = abs(a - b) / max(abs(a), abs(b), colmax[k])
            if math.isnan(rel) or (rel > args.rtol and scaled > args.rtol):
                print(f"FAIL {name}: row {i} (nst={r[0]}) column '{col}': "
                      f"ref={rv} fvode={nv} (scaled diff {scaled:.3e} > {args.rtol:.1e})")
                return False
            max_diff = max(max_diff, scaled)

    if len(rrows) != len(nrows):
        print(f"FAIL {name}: {len(rrows)} reference rows vs {len(nrows)} FVODE rows "
              f"(all {min(len(rrows), len(nrows))} common rows match)")
        return False

    how = "bitwise identical" if max_diff == 0.0 else f"max scaled diff {max_diff:.3e}"
    print(f"PASS {name}: {len(rrows)} rows, same steps and counters, {how}")
    return True


def diverged(name, head, rrows, nrows, real_cols, colmax, row, col, args) -> bool:
    """Step sequences differ (tolerance mode): check the final state instead."""
    rlast, nlast = rrows[-1], nrows[-1]
    ycols = [k for k in real_cols if head[k].startswith("y")]
    worst = max((abs(float(rlast[k]) - float(nlast[k])) / max(colmax[k], 1e-300)
                 for k in ycols), default=0.0)
    nst_r, nst_n = int(rlast[0]), int(nlast[0])
    nst_rel = abs(nst_r - nst_n) / max(nst_r, 1)
    msg = (f"{name}: step sequence diverged at row {row} (column '{col}'); final "
           f"nst {nst_r} vs {nst_n} ({100*nst_rel:.1f}%), final y scaled diff {worst:.2e}")
    if worst <= args.ytol and nst_rel <= args.nsttol:
        print(f"DIVERGED {msg} -- within tolerance (run compare-exact to check fidelity)")
        return True
    print(f"FAIL {msg}")
    return False


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--exact", action="store_true")
    parser.add_argument("--rtol", type=float, default=1e-3)
    parser.add_argument("--ytol", type=float, default=1e-4)
    parser.add_argument("--nsttol", type=float, default=0.10)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    refs = sorted(OUT.glob("trace_*.csv"))
    if not refs:
        print(f"no reference traces in {OUT}; run: pixi run -e reference ref-run")
        return 1

    ok = True
    for ref in refs:
        name = ref.stem.removeprefix("trace_")
        new = OUT / f"fvode_{name}.csv"
        if not new.exists():
            print(f"SKIP {name}: no FVODE trace yet ({new.name})")
            ok = ok and not args.strict
            continue
        ok = compare(name, ref, new, args) and ok
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
