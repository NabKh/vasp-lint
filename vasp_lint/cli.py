"""Command line: vasp-lint [PATH ...]"""

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path

from . import __version__
from .inputs import load
from .rules import ERROR, INFO, RULES, WARNING, run

COLORS = {ERROR: "\033[31m", WARNING: "\033[33m", INFO: "\033[36m"}
RESET = "\033[0m"
ORDER = {ERROR: 0, WARNING: 1, INFO: 2}


def find_calculations(paths):
    """Directories with an INCAR: the path itself, the parent of a file, or every one below."""
    found = []
    for path in map(Path, paths or ["."]):
        if path.is_file():
            path = path.parent
        if (path / "INCAR").is_file():
            found.append(path)
        elif path.is_dir():
            found += sorted(
                p.parent
                for p in path.rglob("INCAR")
                if not any(part.startswith(".") for part in p.relative_to(path).parts)
            )
    unique = []
    for p in found:
        if p.resolve() not in {u.resolve() for u in unique}:
            unique.append(p)
    return unique


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="vasp-lint",
        description="Check VASP input files (INCAR, POSCAR, POTCAR, KPOINTS) for common mistakes before you submit.",
    )
    parser.add_argument("paths", nargs="*", help="calculation folders or files (default: search the current folder)")
    parser.add_argument("--ignore", default="", help="comma-separated rule IDs to skip, e.g. VL001,VL009")
    parser.add_argument("--format", choices=["text", "json"], default="text")
    parser.add_argument("--strict", action="store_true", help="exit with 1 on warnings too")
    parser.add_argument("--no-info", action="store_true", help="hide info-level findings")
    parser.add_argument("--list-rules", action="store_true", help="print all rules and exit")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = parser.parse_args(argv)

    if args.list_rules:
        for r in RULES:
            print(f"{r.id}  {r.severity:7s}  {r.title}\n        {r.url}")
        return 0

    ignore = {s.strip().upper() for s in args.ignore.split(",") if s.strip()}
    calculations = find_calculations(args.paths)
    if not calculations:
        print("vasp-lint: no INCAR found", file=sys.stderr)
        return 2

    report = []
    for folder in calculations:
        calc = load(folder)
        findings = [f for f in run(calc, ignore) if not (args.no_info and f.severity == INFO)]
        findings.sort(key=lambda f: (ORDER[f.severity], f.rule))
        missing = [name for name in ("POSCAR", "POTCAR", "KPOINTS") if getattr(calc, name.lower()) is None]
        report.append((folder, findings, missing))

    if args.format == "json":
        print(
            json.dumps(
                [
                    {"path": str(folder), "missing": missing, "findings": [asdict(f) for f in findings]}
                    for folder, findings, missing in report
                ],
                indent=2,
            )
        )
    else:
        color = sys.stdout.isatty()
        for folder, findings, missing in report:
            for f in findings:
                where = f"{folder / f.file}" + (f":{f.line}" if f.line else "")
                tag = f"{COLORS[f.severity]}{f.severity}{RESET}" if color else f.severity
                print(f"{where}: {f.rule} {tag}: {f.message}\n    {f.url}")
            if missing:
                skipped = "KPOINTS" in missing and "POSCAR" not in missing
                note = " (KSPACING mesh assumed)" if skipped else ""
                print(f"{folder}: note: no readable {', '.join(missing)}; checks that need them were skipped{note}")
        counts = {s: sum(f.severity == s for _, fs, _ in report for f in fs) for s in ORDER}
        print(
            f"\n{len(report)} calculation(s): {counts[ERROR]} error(s), {counts[WARNING]} warning(s), {counts[INFO]} info"
        )

    worst = {f.severity for _, fs, _ in report for f in fs}
    if ERROR in worst or (args.strict and WARNING in worst):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
