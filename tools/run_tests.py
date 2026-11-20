#!/usr/bin/env python3
"""Run the unit tests of every skill's scripts. Standard library only.

Each folder under tests/ is discovered in its own interpreter so that test modules with the
same file name in different skills cannot shadow each other.

usage:
  python3 tools/run_tests.py            # all skills
  python3 tools/run_tests.py NAME ...   # selected skills
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run unit tests for the skill scripts.")
    parser.add_argument("names", nargs="*", help="test folders under tests/ (default: all)")
    parser.add_argument("-q", "--quiet", action="store_true", help="print one line per folder")
    args = parser.parse_args(argv)
    tests = ROOT / "tests"
    folders = [d for d in sorted(tests.iterdir()) if d.is_dir() and not d.name.startswith((".", "_"))]
    if args.names:
        folders = [d for d in folders if d.name in set(args.names)]
    if not folders:
        print("no test folders found", file=sys.stderr)
        return 2
    failed, total = [], 0
    for folder in folders:
        proc = subprocess.run(
            [sys.executable, "-m", "unittest", "discover", "-s", str(folder), "-p", "test_*.py"],
            capture_output=True, text=True, cwd=str(ROOT),
        )
        ran = re.search(r"Ran (\d+) tests? in", proc.stderr)
        count = int(ran.group(1)) if ran else 0
        total += count
        ok = proc.returncode == 0 and count > 0
        print(f"{'ok  ' if ok else 'FAIL'} {folder.name}: {count} tests")
        if not ok:
            failed.append(folder.name)
            if not args.quiet:
                print(proc.stdout[-3000:])
                print(proc.stderr[-6000:])
    print(f"\n{total} tests in {len(folders)} folders; {len(failed)} folder(s) failing" + (f": {', '.join(failed)}" if failed else ""))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
