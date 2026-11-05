"""Compare the output of FVODE example ports with the C SUNDIALS examples.

For each pair below, run the C example (built by ref-build in
reference/build/examples) and the FVODE port (already built by fpm), normalize
both outputs, and diff them:

* lines with workspace lengths are skipped (FVODE does not track workspace),
* case and whitespace are ignored (C prints 1.0e+00, Fortran 1.0E+00),
* Fortran's ".608" is read as "0.608".

Run through pixi (pixi task compare-examples), which builds the FVODE examples
with -ffp-contract=off like the C reference, so that the printed results are
expected to match exactly. Exits non-zero on any difference.
"""

import difflib
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CBIN = ROOT / "reference" / "build" / "examples" / "cvode" / "serial"

# (FVODE example name, C example executable)
PAIRS = [
    ("cv_direct_demo", "cvDirectDemo_ls"),
    ("cv_advdiff_bnd", "cvAdvDiff_bnd"),
]


def normalize(text: str) -> list[str]:
    out = []
    for line in text.splitlines():
        if "workspace length" in line:
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
    for fname, cname in PAIRS:
        c = subprocess.run([str(CBIN / cname)], capture_output=True, text=True, check=True)
        f = subprocess.run(["fpm", "run", "--example", fname, "--profile", "release",
                            "--flag", "-ffp-contract=off"],
                           cwd=ROOT, capture_output=True, text=True, check=True)
        a, b = normalize(c.stdout), normalize(f.stdout)
        diff = [d for d in difflib.unified_diff(a, b, cname, fname, lineterm="", n=0)]
        if diff:
            ok = False
            print(f"FAIL {fname} vs {cname}:")
            print("\n".join(diff[:40]))
        else:
            print(f"PASS {fname}: output identical to {cname} ({len(a)} lines)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
