"""Compare FVODE step traces against the reference C CVODE traces.

For every reference trace  reference/out/trace_<name>.csv  (written by
run_reference.sh), look for the FVODE trace  reference/out/fvode_<name>.csv
(written by the matching FVODE trace program) and compare them row by row:

* integer columns (step counts, orders, counters, return flag) must match exactly;
* real columns must agree to a mixed tolerance (default rtol = 1e-3):
  |a - b| <= rtol * max(|a|, |b|)  or  |a - b| <= rtol * max|column|,
  so values passing through zero are judged against the column's scale.

FVODE is compiled with normal optimizations (e.g. fused multiply-add), so its
real values drift from the C values by rounding, growing slowly over the run.
With --rtol 0 the comparison is exact; use that with an FVODE build made with
-ffp-contract=off (pixi task compare-exact) to locate translation errors.

The first divergent row is reported, which identifies the step where the
translation departs from CVODE. Exits non-zero if any comparison fails.

Usage: python reference/compare.py [--rtol RTOL] [--strict]
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


def compare(name: str, ref: Path, new: Path, rtol: float) -> bool:
    rhead, rrows = read(ref)
    nhead, nrows = read(new)
    if rhead != nhead:
        print(f"FAIL {name}: column headers differ\n  ref: {rhead}\n  new: {nhead}")
        return False

    max_rel = 0.0
    bitwise = True
    colmax = {col: max((abs(float(r[k])) for r in rrows), default=0.0)
              for k, col in enumerate(rhead) if col not in INT_COLS}
    for i, (r, n) in enumerate(zip(rrows, nrows), start=1):
        for col, rv, nv in zip(rhead, r, n):
            if col in INT_COLS:
                if int(rv) != int(nv):
                    print(f"FAIL {name}: row {i} (nst={r[0]}) column '{col}': "
                          f"ref={rv} fvode={nv}")
                    return False
            else:
                a, b = float(rv), float(nv)
                if a != b:
                    bitwise = False
                scale = max(abs(a), abs(b), colmax[col])
                rel = 0.0 if a == b else abs(a - b) / max(abs(a), abs(b))
                mixed = 0.0 if a == b else abs(a - b) / scale
                if math.isnan(rel) or (rel > rtol and mixed > rtol) or (rtol == 0.0 and a != b):
                    print(f"FAIL {name}: row {i} (nst={r[0]}) column '{col}': "
                          f"ref={rv} fvode={nv} (rel diff {rel:.3e} > {rtol:.1e})")
                    return False
                max_rel = max(max_rel, mixed)

    if len(rrows) != len(nrows):
        print(f"FAIL {name}: {len(rrows)} reference steps vs {len(nrows)} FVODE steps "
              f"(all {min(len(rrows), len(nrows))} common rows match)")
        return False

    how = "bitwise identical" if bitwise else f"max scaled diff {max_rel:.3e}"
    print(f"PASS {name}: {len(rrows)} steps, {how}")
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--rtol", type=float, default=1e-3)
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
        ok = compare(name, ref, new, args.rtol) and ok
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
