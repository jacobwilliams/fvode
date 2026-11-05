"""Compare the output of FVODE example ports with the C SUNDIALS examples.

For each pair below, run the C example (built by ref-build in
reference/build/examples) and the FVODE port (already built by fpm), normalize
both outputs, and diff them:

* lines with workspace lengths (or lenrw/leniw) are skipped (FVODE does not track workspace),
* case and whitespace are ignored (C prints 1.0e+00, Fortran 1.0E+00),
* Fortran's ".608" is read as "0.608".

Run through pixi (pixi task compare-examples), which builds the FVODE examples
with -ffp-contract=off like the C reference, so that the printed results are
expected to match exactly. Exits non-zero on any difference.

A pair may allow a relative tolerance `rtol` on real numbers (integers, i.e.
all step counts and statistics, must always match exactly). This is needed
where FVODE's integer powers (x**i, by repeated squaring) round differently
from SUNDIALS's SUNRpowerI (a multiplication loop): the steps are the same,
but the printed values can differ in their last digits.
"""

import difflib
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CBIN = ROOT / "reference" / "build" / "examples" / "cvode" / "serial"

# (FVODE example name, C example executable, relative tolerance on reals[,
# command-line arguments for both programs])
PAIRS = [
    ("cv_direct_demo", "cvDirectDemo_ls", 0.0),
    ("cv_advdiff_bnd", "cvAdvDiff_bnd", 0.0),
    ("cv_roberts_dns", "cvRoberts_dns", 0.0),
    ("cv_roberts_dns_constraints", "cvRoberts_dns_constraints", 0.0),
    # x**i vs SUNRpowerI: last two outputs differ by ~1e-12 absolute (identical
    # when FVODE uses a multiplication loop); all statistics identical
    ("cv_roberts_dns_negsol", "cvRoberts_dns_negsol", 1.0e-4),
    ("cv_rocket_dns", "cvRocket_dns", 0.0),
    ("cv_disc_dns", "cvDisc_dns", 0.0),
    ("cv_diurnal_kry", "cvDiurnal_kry", 0.0),
    ("cv_diurnal_kry_bp", "cvDiurnal_kry_bp", 0.0),
    ("cv_krylov_demo_ls", "cvKrylovDemo_ls", 0.0),
    ("cv_krylov_demo_ls", "cvKrylovDemo_ls", 0.0, ["1"]),  # nrmfac = sqrt(N) given
    ("cv_krylov_demo_ls", "cvKrylovDemo_ls", 0.0, ["2"]),  # nrmfac from a dot product
    ("cv_krylov_demo_prec", "cvKrylovDemo_prec", 0.0),
]

NUMBER = re.compile(r"^[-+]?(\d+\.\d*|\.\d+|\d+)(e[-+]?\d+)?$")


def same_line(a: str, b: str, rtol: float) -> bool:
    """Exact match, or equal up to rtol in the real-valued tokens."""
    if a == b:
        return True
    ta, tb = a.split(" "), b.split(" ")
    if len(ta) != len(tb) or rtol == 0.0:
        return False
    for x, y in zip(ta, tb):
        if x == y:
            continue
        if not (NUMBER.match(x) and NUMBER.match(y)) or not any(c in x + y for c in ".e"):
            return False  # non-numeric or integer tokens must match exactly
        u, v = float(x), float(y)
        if abs(u - v) > rtol * max(abs(u), abs(v)):
            return False
    return True


def normalize(text: str) -> list[str]:
    out = []
    for line in text.splitlines():
        if "workspace length" in line or re.match(r"\s*leniw|\s*lenrw", line):
            continue
        line = line.strip().lower()
        if not line:
            continue
        # fpm's own build messages
        if re.match(r"^\[\s*\d+%\]", line) or line.startswith("project "):
            continue
        line = re.sub(r"(?<![\d.])\.(\d)", r"0.\1", line)
        out.append(re.sub(r"\s+", " ", line))
    return out


def main() -> int:
    ok = True
    for pair in PAIRS:
        fname, cname, rtol = pair[:3]
        args = list(pair[3]) if len(pair) > 3 else []
        c = subprocess.run([str(CBIN / cname), *args], cwd=ROOT / "reference" / "out",
                           capture_output=True, text=True, check=True)
        f = subprocess.run(["fpm", "run", "--example", fname, "--profile", "release",
                            "--flag", "-ffp-contract=off", *(["--", *args] if args else [])],
                           cwd=ROOT, capture_output=True, text=True, check=True)
        if args:
            fname, cname = f"{fname} {' '.join(args)}", f"{cname} {' '.join(args)}"
        a, b = normalize(c.stdout), normalize(f.stdout)
        if len(a) == len(b) and all(same_line(x, y, rtol) for x, y in zip(a, b)):
            ndiff = sum(x != y for x, y in zip(a, b))
            if ndiff == 0:
                print(f"PASS {fname}: output identical to {cname} ({len(a)} lines)")
            else:
                print(f"PASS {fname}: output matches {cname} ({len(a)} lines, "
                      f"{ndiff} within rtol = {rtol:g})")
        else:
            ok = False
            diff = [d for d in difflib.unified_diff(a, b, cname, fname, lineterm="", n=0)]
            print(f"FAIL {fname} vs {cname}:")
            print("\n".join(diff[:40]))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
