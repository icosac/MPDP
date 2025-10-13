from __future__ import annotations

import argparse
import ast
import cProfile
import importlib
import io
import pstats
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence


PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parent


def _load_callable(target: str) -> Callable[..., Any]:
    module_path, _, attr = target.rpartition(".")
    if not module_path:
        raise ValueError(f"Target '{target}' must be a dotted path to a callable")
    module = importlib.import_module(module_path)
    try:
        candidate = getattr(module, attr)
    except AttributeError as exc:  # pragma: no cover - defensive
        raise ValueError(f"Module '{module_path}' has no attribute '{attr}'") from exc
    if not callable(candidate):
        raise TypeError(f"Target '{target}' is not callable")
    return candidate


def _format_location(filename: str, lineno: int, func_name: str) -> str:
    if filename.startswith("<"):
        return f"{filename}:{lineno}:{func_name}"
    path = Path(filename)
    try:
        resolved = path.resolve()
    except (OSError, RuntimeError, ValueError):
        resolved = path
    try:
        rel_path = resolved.relative_to(REPO_ROOT.resolve())
        display = rel_path.as_posix()
    except ValueError:
        display = resolved.name
    return f"{display}:{lineno}:{func_name}"


def _iter_profile_rows(stats: pstats.Stats, include_external: bool) -> Iterable[tuple[str, int, int, float, float, float]]:
    package_marker = str(PACKAGE_ROOT.resolve())
    func_list = getattr(stats, "fcn_list", list(stats.stats.keys()))
    for func_identifier in func_list:
        filename, lineno, func_name = func_identifier
        try:
            resolved = Path(filename).resolve()
        except (OSError, RuntimeError, ValueError):
            resolved = Path(filename)
        if not include_external and package_marker not in str(resolved):
            continue
        ccalls, ncalls, tt, ct, _ = stats.stats[func_identifier]
        per_call = ct / ncalls if ncalls else 0.0
        location = _format_location(filename, lineno, func_name)
        yield location, ccalls, ncalls, tt, ct, per_call


def render_stats_table(stats: pstats.Stats, *, include_external: bool = False, limit: int | None = 30) -> str:
    header = f"{'ncalls':>12} {'prim':>6} {'tottime':>10} {'cumtime':>10} {'percall':>10}  location"
    separator = "-" * len(header)
    lines = [header, separator]
    count = 0
    for location, ccalls, ncalls, tt, ct, per_call in _iter_profile_rows(stats, include_external):
        ncalls_display = f"{ncalls}" if ccalls == ncalls else f"{ncalls}/{ccalls}"
        lines.append(f"{ncalls_display:>12} {ccalls:>6} {tt:>10.6f} {ct:>10.6f} {per_call:>10.6f}  {location}")
        count += 1
        if limit is not None and count >= limit:
            break
    if count == 0:
        lines.append("(no entries matched the current filter)")
    return "\n".join(lines)


def profile_callable(
    func: Callable[..., Any],
    *args: Any,
    sort_by: str = "cumulative",
    limit: int | None = 30,
    include_external: bool = False,
    strip_dirs: bool = True,
    **kwargs: Any,
) -> tuple[Any, str, pstats.Stats]:
    profiler = cProfile.Profile()
    profiler.enable()
    try:
        result = func(*args, **kwargs)
    finally:
        profiler.disable()
    stats = pstats.Stats(profiler, stream=io.StringIO())
    if strip_dirs:
        stats.strip_dirs()
    stats.sort_stats(sort_by)
    table = render_stats_table(stats, include_external=include_external, limit=limit)
    return result, table, stats


def _parse_positional(value: str) -> tuple[Any, ...]:
    parsed = ast.literal_eval(value)
    if isinstance(parsed, tuple):
        return parsed
    if isinstance(parsed, list):
        return tuple(parsed)
    raise ValueError(f"Expected list or tuple literal for positional arguments, got {type(parsed).__name__}")


def _parse_keyword(value: str) -> dict[str, Any]:
    parsed = ast.literal_eval(value)
    if isinstance(parsed, dict):
        return parsed
    raise ValueError(f"Expected dict literal for keyword arguments, got {type(parsed).__name__}")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Profile a callable within pympdp using cProfile.")
    parser.add_argument("target", help="Dotted path to callable, e.g. pympdp.dp.dp.DP.solve_dp")
    parser.add_argument("--args", default="[]", help="Literal list/tuple of positional arguments (default: [])")
    parser.add_argument("--kwargs", default="{}", help="Literal dict of keyword arguments (default: {})")
    parser.add_argument("--sort", default="cumulative", help="Sort key for pstats (default: cumulative)")
    parser.add_argument("--limit", type=int, default=30, help="Maximum number of rows in the table (default: 30)")
    parser.add_argument("--no-strip-dirs", action="store_true", help="Disable directory stripping in pstats output")
    parser.add_argument("--include-external", action="store_true", help="Include non-pympdp entries in the table")
    args = parser.parse_args(argv)

    callable_obj = _load_callable(args.target)
    positional = _parse_positional(args.args)
    keyword = _parse_keyword(args.kwargs)

    _, table, _ = profile_callable(
        callable_obj,
        *positional,
        sort_by=args.sort,
        limit=args.limit,
        include_external=args.include_external,
        strip_dirs=not args.no_strip_dirs,
        **keyword,
    )
    print(table)
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI utility
    raise SystemExit(main())
