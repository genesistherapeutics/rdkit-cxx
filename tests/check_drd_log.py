# Copyright (c) 2026 Steven Kearnes
# All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

"""Fails on races in RDKit or the wrappers in a log from Valgrind's DRD.

CI runs tests/threads.rs under DRD (CONVENTIONS.md U4). DRD also reports Rust's
own thread bookkeeping, which synchronizes with atomics it cannot see, so this
script fails only on reports that involve RDKit: an access to an RDKit
library's static data, or a conflicting access made from RDKit or wrapper code.
"""

import argparse
import re
import sys

PREFIX = re.compile(r"^==\d+== ?")
FRAME = re.compile(r"^\s*(at|by) 0x[0-9A-F]+: ")
RDKIT_DATA = re.compile(r"(BSS|Data) section of \S*/libRDKit")
RDKIT_CODE = re.compile(r"RDKit::|rdkit_cxx::|\(in \S*/libRDKit")


def reports(lines: list[str]) -> list[list[str]]:
    """Splits the log at blank lines into DRD's error reports.

    Each report starts with its description, after a "Thread N:" line when the
    reporting thread changes, and is followed by the stack that made it.
    """
    paragraphs = [[]]
    for line in lines:
        if line:
            paragraphs[-1].append(line)
        elif paragraphs[-1]:
            paragraphs.append([])
    return [p for p in paragraphs if any(FRAME.match(line) for line in p)]


def involves_rdkit(report: list[str]) -> bool:
    if any(RDKIT_DATA.search(line) for line in report):
        return True
    # The stack of the access DRD reports (not those of the other thread's
    # segment) is the run of frames right after the description.
    start = 2 if re.match(r"Thread \d+:$", report[0]) else 1
    access = []
    for line in report[start:]:
        if not FRAME.match(line):
            break
        access.append(line)
    return any(RDKIT_CODE.search(line) for line in access)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", help="the file DRD wrote with --log-file")
    args = parser.parse_args()
    with open(args.log) as f:
        lines = [PREFIX.sub("", line.rstrip("\n")) for line in f]
    if not any(line.startswith("ERROR SUMMARY") for line in lines):
        print(f"{args.log} is not a complete DRD log", file=sys.stderr)
        return 1
    races = [report for report in reports(lines) if involves_rdkit(report)]
    for report in races:
        print("\n".join(report))
    print(f"{len(races)} DRD reports involve RDKit or the wrappers", file=sys.stderr)
    return 1 if races else 0


if __name__ == "__main__":
    sys.exit(main())
