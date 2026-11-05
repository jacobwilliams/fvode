"""Check that the Fortran sources are documented (FORD style).

Rules:
1. every procedure (subroutine/function, including interface bodies and
   module procedures) has a `!>` comment block before it or a `!!` docstring
   right after its first statement;
2. every derived type has a `!>` block before it or a `!!` docstring after
   its `type ::` line;
3. every declaration (arguments, local variables, parameters, module
   variables, type components, procedure pointers, type-bound procedures)
   has a `!!` comment, on the same line or on the following line(s).

Usage: python tools/doc_audit.py [paths...]   (default: src example test reference/fortran)
Prints one line per finding and exits non-zero if there are any.
"""

import re
import sys
from pathlib import Path

DECL = re.compile(
    r"^\s*(integer|real|logical|character|complex|type\s*\(|class\s*\(|procedure\s*(\(|,|::|\s))"
    r"[^!]*::", re.I)
PROC = re.compile(
    r"^\s*((pure|elemental|impure|recursive|module|non_recursive)\s+)*"
    r"((integer|real|logical|character)(\s*\([^)]*\))?\s+)?"
    r"(subroutine|function)\s+\w+", re.I)
TYPEDEF = re.compile(r"^\s*type\s*(,[^:]*)?::\s*\w+", re.I)
END = re.compile(r"^\s*end\b", re.I)


def logical_lines(lines):
    """Join continuation lines; yield (first line number, text, next index)."""
    i = 0
    while i < len(lines):
        start = i
        text = lines[i].rstrip("\n")
        code = text.split("!")[0] if "!!" not in text else text.split("!!")[0]
        while code.rstrip().endswith("&") and i + 1 < len(lines):
            i += 1
            text += "\n" + lines[i].rstrip("\n")
            nxt = lines[i]
            code = nxt.split("!!")[0]
        yield start, text, i + 1
        i += 1


def has_doc_after(lines, j):
    """Is line j (after a statement) a `!!` docstring line?"""
    return j < len(lines) and lines[j].lstrip().startswith("!!")


def has_doc_before(lines, i):
    """Is there a `!>` block above line i (skipping separators and blank lines)?"""
    k = i - 1
    while k >= 0:
        s = lines[k].strip()
        if s.startswith("!>"):
            return True
        if s == "" or re.fullmatch(r"!\*+|!=+|!-+", s):
            k -= 1
            continue
        if s.startswith("!"):     # part of a !> block (continuation lines)
            k -= 1
            continue
        return False
    return False


def audit(path: Path) -> list[str]:
    lines = path.read_text().splitlines()
    found = []
    for start, text, nxt in logical_lines(lines):
        first = text.split("\n")[0]
        stripped = first.strip()
        if stripped.startswith("!") or not stripped:
            continue
        if END.match(first):
            continue
        documented_inline = "!!" in text
        if PROC.match(first) and not re.match(r"^\s*procedure\b", first, re.I):
            if not (has_doc_before(lines, start) or has_doc_after(lines, nxt)):
                found.append(f"{path}:{start+1}: procedure without docstring: {stripped[:70]}")
        elif TYPEDEF.match(first):
            if not (has_doc_before(lines, start) or has_doc_after(lines, nxt) or documented_inline):
                found.append(f"{path}:{start+1}: type without docstring: {stripped[:70]}")
        elif DECL.match(first):
            if not (documented_inline or has_doc_after(lines, nxt)):
                found.append(f"{path}:{start+1}: declaration without !! comment: {stripped[:70]}")
    return found


def main() -> int:
    roots = [Path(p) for p in (sys.argv[1:] or ["src", "example", "test", "reference/fortran"])]
    files = sorted(f for r in roots for f in ([r] if r.is_file() else r.rglob("*.[fF]90")))
    findings = [x for f in files for x in audit(f)]
    for x in findings:
        print(x)
    print(f"{len(findings)} finding(s) in {len(files)} file(s)")
    return 1 if findings else 0


if __name__ == "__main__":
    sys.exit(main())
