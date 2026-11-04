#!/usr/bin/env python3
"""Compare results collected on different boards.

    python3 compare_boards.py OUR_FOLDER THEIR_FOLDER [MORE...]
    python3 compare_boards.py UNO-Q=. beagle="~/Downloads/bms_result_beaglebone"

Each folder is either a `collect.sh` output (run_1.txt ... run_3.txt, data/,
environment.txt) or the folder this script lives in (run_output.txt, data/).
Prefix a folder with `LABEL=` to name it by hand.

The point of the exercise is that the physics must come out identical and only
the timing may differ. This script enforces that: it fails loudly if any board
reports a different trip step, temperature or energy.
"""
import re
import sys
from pathlib import Path

FIELDS = {
    "step": (r"== LIMIT REACHED\s+t = [\d.]+ s\s+\(step (\d+)\)", int),
    "t": (r"== LIMIT REACHED\s+t = ([\d.]+) s", float),
    "surface": (r"surface\s+:\s+([\d.]+) C", float),
    "core": (r"centre\s+:\s+([\d.]+) C", float),
    "gap": (r"\(([+-][\d.]+) K above the surface\)", float),
    "stored": (r"Energy stored above amb : ([\d.]+) J", float),
    "generated": (r"Energy generated\s+: ([\d.]+) J", float),
    "convected": (r"Energy convected away\s+: ([\d.]+) J", float),
    "wall": (r"wall\s+:\s+([\d.]+) s", float),
    "per_step": (r"per step:\s+([\d.]+) ms", float),
    "workload": (r"workload: (.*)", str),
}

PHYSICS = ("step", "t", "surface", "core", "gap", "stored", "generated", "convected")


def parse(path):
    text = path.read_text(errors="replace")
    out = {}
    for name, (pattern, cast) in FIELDS.items():
        m = re.search(pattern, text)
        if m:
            out[name] = cast(m.group(1).strip())
    return out


def runs_in(folder):
    numbered = sorted(folder.glob("run_[123].txt"))
    if numbered:
        return numbered
    single = folder / "run_output.txt"
    return [single] if single.exists() else []


def label(folder):
    env = folder / "environment.txt"
    if env.exists():
        text = env.read_text(errors="replace")
        # lscpu says "Model name:", /proc/cpuinfo says "model name", ARM boards
        # often only have "Hardware". Skip a bare "Model:" (that one is a number).
        for pattern in (r"(?im)^\s*Model name\s*:\s*(.+)$",
                        r"(?im)^\s*Hardware\s*:\s*(.+)$",
                        r"(?im)^\s*Processor\s*:\s*(.+)$"):
            m = re.search(pattern, text)
            if m and m.group(1).strip():
                return "%s (%s)" % (folder.name, m.group(1).strip()[:40])
    return folder.name


def main(argv):
    if len(argv) < 3:
        print(__doc__)
        return 2

    boards = []
    for arg in argv[1:]:
        # "LABEL=PATH" lets a folder without environment.txt be named by hand,
        # e.g. "UNO-Q=." ; a bare path is labelled from environment.txt.
        if "=" in arg:
            forced, _, arg = arg.partition("=")
            forced = forced.strip()
        else:
            forced = None
        folder = Path(arg).resolve()
        files = runs_in(folder)
        if not files:
            print("no run_*.txt or run_output.txt in %s" % folder)
            return 2
        parsed = [parse(f) for f in files]
        plain = [p for p in parsed if p.get("workload", "").startswith("plain")]
        if not plain:
            plain = parsed
        timings = [p["wall"] for p in plain if "wall" in p]
        boards.append({
            "name": forced or label(folder),
            "folder": folder,
            "runs": parsed,
            "plain": plain[0],
            "timings": timings,
        })

    # ---- the physics must agree -------------------------------------------
    ref = boards[0]["plain"]
    problems = []
    for b in boards:
        for p in b["runs"]:
            for key in PHYSICS:
                if key not in p or key not in ref:
                    continue
                a, r = p[key], ref[key]
                if key == "step":
                    if a != r:
                        problems.append("%s: step %s != %s" % (b["name"], a, r))
                elif abs(a - r) > 1e-6:
                    problems.append("%s: %s %s != %s" % (b["name"], key, a, r))

    print("=" * 78)
    print("physics")
    print("=" * 78)
    print("  trip step      %d" % ref["step"])
    print("  trip time      %.4f s" % ref["t"])
    print("  surface        %.4f C" % ref["surface"])
    print("  core           %.4f C" % ref["core"])
    print("  core-surface   %+.4f K" % ref["gap"])
    print("  stored         %.4f J" % ref["stored"])
    if problems:
        print()
        print("  MISMATCH — these boards are not running the same model:")
        for p in problems:
            print("    " + p)
        print("  Timing comparison below is NOT meaningful.")
    else:
        print("  identical on all %d boards" % len(boards))

    # ---- timing ------------------------------------------------------------
    print()
    print("=" * 78)
    print("timing (plain run, no export, no stream)")
    print("=" * 78)
    header = "  %-42s %10s %10s %14s" % ("board", "mean s", "ms/step", "vs first")
    print(header)
    print("  " + "-" * (len(header) - 2))
    base = None
    for b in boards:
        if not b["timings"]:
            print("  %-42s %10s" % (b["name"][:42], "n/a"))
            continue
        mean = sum(b["timings"]) / len(b["timings"])
        # Derived from the mean, not from one run's own printout, so the two
        # columns describe the same set of runs.
        steps = b["plain"].get("step")
        per_step = mean / steps * 1000.0 if steps else float("nan")
        if base is None:
            base = mean
            ratio = "1.00x"
        elif 0.97 <= mean / base <= 1.03:
            # Three runs on one board already vary by more than this.
            ratio = "~same"
        elif mean < base:
            ratio = "%.1fx faster" % (base / mean)
        else:
            ratio = "%.2fx slower" % (mean / base)
        runs = " ".join("%.3f" % t for t in b["timings"])
        print("  %-42s %10.4f %10.4f %14s   [%s]"
              % (b["name"][:42], mean, per_step, ratio, runs))

    # ---- exported data, if both sides ran --export -------------------------
    with_data = [b for b in boards if (b["folder"] / "data" / "temperature_history.csv").exists()]
    if len(with_data) >= 2:
        print()
        print("=" * 78)
        print("exported field")
        print("=" * 78)
        a = with_data[0]["folder"] / "data"
        for b in with_data[1:]:
            d = b["folder"] / "data"
            same_bytes = (a / "midplane_slices.csv").read_bytes() == (d / "midplane_slices.csv").read_bytes()
            print("  %-42s slices %s" % (b["name"][:42],
                                         "byte-identical" if same_bytes else "differ (check precision)"))

    print()
    print("Only the timing is supposed to differ. Anything else means the")
    print("boards are not running the same model.")
    return 0 if not problems else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
