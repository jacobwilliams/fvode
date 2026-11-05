"""Compare FVODE step traces against the reference C CVODE traces.

For every reference trace  reference/out/trace_<name>.csv  (written by
run_reference.sh), look for the FVODE trace  reference/out/fvode_<name>.csv
(written by the matching FVODE trace program) and compare them row by row:

* integer columns (step counts, orders, counters, return flag) must match exactly;
* real columns must agree to a relative tolerance (default 1e-12, see --rtol).

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
                scale = max(abs(a), abs(b))
                rel = 0.0 if scale == 0.0 else abs(a - b) / scale
                if math.isnan(rel) or rel > rtol:
                    print(f"FAIL {name}: row {i} (nst={r[0]}) column '{col}': "
                          f"ref={rv} fvode={nv} (rel diff {rel:.3e} > {rtol:.1e})")
                    return False
                max_rel = max(max_rel, rel)

    if len(rrows) != len(nrows):
        print(f"FAIL {name}: {len(rrows)} reference steps vs {len(nrows)} FVODE steps "
              f"(all {min(len(rrows), len(nrows))} common rows match)")
        return False

    how = "bitwise identical" if bitwise else f"max rel diff {max_rel:.3e}"
    print(f"PASS {name}: {len(rrows)} steps, {how}")
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--rtol", type=float, default=1e-12)
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
