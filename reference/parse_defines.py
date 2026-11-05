"""Print the numeric #define constants local to CVODE .c files as "NAME VALUE".

Usage: python -I parse_defines.py <sundials-source-dir>

Header constants are printed by constants.c. These constants live only in .c
files, so they cannot be #included and are parsed from the source text.
"""
import re
import sys
from pathlib import Path

FILES = ["cvode.c", "cvode_nls.c", "cvode_ls.c", "cvode_diag.c", "cvode_bandpre.c"]
DEFINE = re.compile(r"^#define\s+([A-Z][A-Z0-9_]*)\s+(.+?)\s*(?:/\*.*)?$")
RCONST = re.compile(r"^SUN_RCONST\(\s*([-+0-9.eE]+)\s*\)$")
INT = re.compile(r"^[-+]?\d+$")


def main(src: Path) -> None:
    values: dict[str, str] = {}
    for name in FILES:
        for line in (src / "src" / "cvode" / name).read_text().splitlines():
            m = DEFINE.match(line)
            if not m:
                continue
            key, val = m.group(1), m.group(2).strip()
            if (r := RCONST.match(val)):
                out = repr(float(r.group(1)))
            elif INT.match(val):
                out = str(int(val))
            elif val in values:          # alias, e.g. H_BIAS HALF
                out = values[val]
            else:
                continue                 # message strings, macros, etc.
            if key in values and values[key] != out:
                sys.exit(f"conflicting values for {key}: {values[key]} vs {out}")
            values[key] = out
    for key, out in values.items():
        print(key, out)


if __name__ == "__main__":
    main(Path(sys.argv[1]))
