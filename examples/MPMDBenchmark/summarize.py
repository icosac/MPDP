#!/usr/bin/env python3
"""Merge the per-solver result files into one comparison report.

Reads the problem file and every ``results_*.csv`` produced by the drivers and
writes a Markdown summary plus a merged CSV. Standard library only, so it runs
with whatever python3 is on the machine.

Usage:
    summarize.py --problems <path> --results <dir> --profile <name> [--out <dir>]
"""

import argparse
import csv
import math
import os
import sys

#: Solver used as the reference for correctness and for the speed-up column.
BASELINE = "cpu"


def read_problems(path):
    """Return ``{id: problem dict}`` from the shared problem file."""
    problems = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            problems[int(row["id"])] = {
                "tag": row["tag"],
                "n_points": int(row["n_points"]),
                "kmax": float(row["kmax"]),
                "discr": int(row["discr"]),
                "nref": int(row["nref"]),
            }
    return problems


def read_results(path):
    """Return ``{id: result dict}`` from one solver's result file."""
    results = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            angles = [float(a) for a in row["angles"].split(";") if a]
            results[int(row["id"])] = {
                "solver": row["solver"],
                "status": row["status"],
                "length": float(row["length"]),
                "length_check": float(row["length_check"]),
                "time_ms": float(row["time_ms"]),
                "time_solver_ms": float(row["time_solver_ms"]),
                "angles": angles,
            }
    return results


def angle_gap(a, b):
    """Largest difference between two angle lists, modulo a full turn."""
    if len(a) != len(b) or not a:
        return float("nan")
    worst = 0.0
    for x, y in zip(a, b):
        d = math.fmod(abs(x - y), 2.0 * math.pi)
        d = min(d, 2.0 * math.pi - d)
        worst = max(worst, d)
    return worst


def fmt(value, spec="{:.3g}"):
    """Format a float, tolerating NaN."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "n/a"
    return spec.format(value)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--problems", required=True)
    parser.add_argument("--results", required=True, help="directory holding results_*.csv")
    parser.add_argument("--profile", required=True)
    parser.add_argument("--out", default=None, help="defaults to --results")
    args = parser.parse_args()

    out_dir = args.out or args.results
    problems = read_problems(args.problems)

    prefix = "results_{}_".format(args.profile)
    files = sorted(
        f for f in os.listdir(args.results)
        if f.startswith(prefix) and f.endswith(".csv")
    )
    if not files:
        sys.exit("no result files matching {}*.csv in {}".format(prefix, args.results))

    solvers = {}
    for name in files:
        rows = read_results(os.path.join(args.results, name))
        if not rows:
            continue
        label = next(iter(rows.values()))["solver"]
        solvers[label] = rows

    if BASELINE not in solvers:
        sys.exit("baseline solver '{}' not found among {}".format(
            BASELINE, ", ".join(sorted(solvers))))

    base = solvers[BASELINE]
    others = [s for s in sorted(solvers) if s != BASELINE]

    merged_path = os.path.join(out_dir, "merged_{}.csv".format(args.profile))
    with open(merged_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "id", "tag", "n_points", "discr", "nref", "solver", "status",
            "length", "length_check", "self_gap", "abs_err_vs_cpu",
            "rel_err_vs_cpu", "max_angle_gap", "time_ms", "time_solver_ms",
            "speedup_vs_cpu",
        ])
        for pid in sorted(problems):
            prob = problems[pid]
            for label in [BASELINE] + others:
                res = solvers[label].get(pid)
                if res is None:
                    continue
                ref = base.get(pid)
                ok = res["status"] == "ok" and ref is not None and ref["status"] == "ok"
                abs_err = res["length_check"] - ref["length_check"] if ok else float("nan")
                rel_err = abs_err / ref["length_check"] if ok and ref["length_check"] else float("nan")
                gap = angle_gap(res["angles"], ref["angles"]) if ok else float("nan")
                speed = ref["time_ms"] / res["time_ms"] if ok and res["time_ms"] > 0 else float("nan")
                self_gap = abs(res["length"] - res["length_check"]) if res["status"] == "ok" else float("nan")
                writer.writerow([
                    pid, prob["tag"], prob["n_points"], prob["discr"], prob["nref"],
                    label, res["status"],
                    repr(res["length"]), repr(res["length_check"]), repr(self_gap),
                    repr(abs_err), repr(rel_err), repr(gap),
                    repr(res["time_ms"]), repr(res["time_solver_ms"]), repr(speed),
                ])

    lines = []
    lines.append("# MPMDBenchmark comparison - profile `{}`".format(args.profile))
    lines.append("")
    lines.append("Baseline: `{}` (srcCC). `length` is what the solver's own DP reports, "
                 "`check` is the path length recomputed from the returned angles.".format(BASELINE))
    lines.append("")

    lines.append("## Per instance")
    lines.append("")
    header = ("| id | tag | N | discr | solver | status | length | check | "
              "self gap | rel. err vs cpu | max angle gap | time (ms) | speed-up |")
    lines.append(header)
    lines.append("|" + "---|" * 13)
    for pid in sorted(problems):
        prob = problems[pid]
        for label in [BASELINE] + others:
            res = solvers[label].get(pid)
            if res is None:
                continue
            ref = base.get(pid)
            ok = res["status"] == "ok" and ref is not None and ref["status"] == "ok"
            abs_err = res["length_check"] - ref["length_check"] if ok else float("nan")
            rel_err = abs_err / ref["length_check"] if ok and ref["length_check"] else float("nan")
            gap = angle_gap(res["angles"], ref["angles"]) if ok else float("nan")
            speed = ref["time_ms"] / res["time_ms"] if ok and res["time_ms"] > 0 else float("nan")
            self_gap = abs(res["length"] - res["length_check"]) if res["status"] == "ok" else float("nan")
            lines.append("| {} | {} | {} | {} | `{}` | {} | {} | {} | {} | {} | {} | {} | {} |".format(
                pid, prob["tag"], prob["n_points"], prob["discr"], label, res["status"],
                fmt(res["length"], "{:.10g}"), fmt(res["length_check"], "{:.10g}"),
                fmt(self_gap, "{:.2e}"), fmt(rel_err, "{:+.2e}"), fmt(gap, "{:.2e}"),
                fmt(res["time_ms"], "{:.2f}"),
                fmt(speed, "{:.1f}x") if not math.isnan(speed) else "n/a"))
    lines.append("")

    lines.append("## Totals")
    lines.append("")
    lines.append("| solver | solved | skipped | failed | total time (s) | speed-up vs cpu | "
                 "worst rel. err | worst self gap |")
    lines.append("|" + "---|" * 8)
    for label in [BASELINE] + others:
        rows = solvers[label]
        solved = [r for r in rows.values() if r["status"] == "ok"]
        skipped = [r for r in rows.values() if r["status"] == "skipped"]
        failed = [r for r in rows.values() if r["status"].startswith("error")]
        total = sum(r["time_ms"] for r in solved) / 1000.0
        common = [pid for pid, r in rows.items()
                  if r["status"] == "ok" and base.get(pid, {}).get("status") == "ok"]
        base_total = sum(base[pid]["time_ms"] for pid in common) / 1000.0
        mine_total = sum(rows[pid]["time_ms"] for pid in common) / 1000.0
        speed = base_total / mine_total if mine_total > 0 else float("nan")
        errs = [abs(rows[pid]["length_check"] - base[pid]["length_check"]) /
                base[pid]["length_check"]
                for pid in common if base[pid]["length_check"]]
        gaps = [abs(rows[pid]["length"] - rows[pid]["length_check"]) for pid in common]
        lines.append("| `{}` | {} | {} | {} | {} | {} | {} | {} |".format(
            label, len(solved), len(skipped), len(failed), fmt(total, "{:.3f}"),
            fmt(speed, "{:.1f}x"), fmt(max(errs) if errs else float("nan"), "{:.2e}"),
            fmt(max(gaps) if gaps else float("nan"), "{:.2e}")))
    lines.append("")
    lines.append("Merged per-row data: `{}`".format(os.path.basename(merged_path)))
    lines.append("")

    report = "\n".join(lines)
    report_path = os.path.join(out_dir, "summary_{}.md".format(args.profile))
    with open(report_path, "w") as handle:
        handle.write(report)

    print(report)
    print("written: {}".format(report_path))
    print("written: {}".format(merged_path))


if __name__ == "__main__":
    main()
