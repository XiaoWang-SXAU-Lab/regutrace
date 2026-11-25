#!/usr/bin/env python3
"""Run a command and record its peak resident memory and wall time.

Windows uses `tasklist`; POSIX uses /proc/<pid>/status VmHWM. The wrapper is
deliberately external so the measured binary does not pay instrumentation
cost.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path


def windows_memory_kb(pid: int) -> int | None:
    try:
        out = subprocess.check_output(
            ["tasklist", "/FI", f"PID eq {pid}", "/FO", "CSV", "/NH"],
            stderr=subprocess.DEVNULL,
            text=True,
            timeout=5,
        )
    except Exception:
        return None
    line = out.strip().splitlines()
    if not line or line[0].startswith("INFO:"):
        return None
    fields = next(iter(__import__("csv").reader(line, delimiter=",")), [])
    if len(fields) < 5:
        return None
    value = fields[4].replace(",", "").replace("K", "").strip()
    try:
        return int(value)
    except ValueError:
        return None


def posix_memory_kb(pid: int) -> int | None:
    try:
        for line in Path(f"/proc/{pid}/status").read_text().splitlines():
            if line.startswith("VmHWM:"):
                return int(line.split()[1])
    except OSError:
        return None
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--log", required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        parser.error("a command is required after --")

    started = time.time()
    with open(args.log, "w", encoding="utf-8") as log:
        process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
        peak = 0
        while process.poll() is None:
            memory = (
                windows_memory_kb(process.pid)
                if os.name == "nt"
                else posix_memory_kb(process.pid)
            )
            if memory is not None:
                peak = max(peak, memory)
            time.sleep(0.25)
        code = process.returncode
    elapsed = time.time() - started
    result = {
        "label": args.label,
        "exit_code": code,
        "wall_seconds": elapsed,
        "peak_working_set_kb": peak,
        "peak_working_set_mib": peak / 1024.0,
        "command": command,
    }
    Path(args.out).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0 if code == 0 else code


if __name__ == "__main__":
    sys.exit(main())
