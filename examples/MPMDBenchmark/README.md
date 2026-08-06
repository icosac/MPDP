# MPMDBenchmark

A comparison harness for the multi-point Dubins solvers: the CPU reference in
`srcCC` (`DP::solveDP`) against the GPU solver in `srcCU`
(`mpdp::gpu::solveDP`).

It was originally built to compare three implementations - the CPU, the stale
GPU code that used to live in `srcCU`, and its replacement. The old GPU solver
has since been removed, so what remains is a CPU-vs-GPU correctness and
throughput check.

Everything it produces goes into `examples/MPMDBenchmark/results/`. The only thing written
outside this directory is the CMake build tree, which defaults to
`<repo>/build-mpmd-benchmark`.

## Running it

```bash
./examples/MPMDBenchmark/run_benchmark.sh
```

That configures, builds, generates the problem set, runs every solver on it and
writes the report. Useful options:

| option | default | meaning |
|---|---|---|
| `--budget <n>` | `1000000` | total point-to-point curves the `budget` profile aims for |
| `--max-seconds <s>` | `60` | wall-clock cap per solver per profile |
| `--instances <n>` | `1` | random instances per configuration |
| `--reps <n>` | `3` | timing repetitions per instance; the fastest is reported |
| `--profiles <list>` | `budget,stress` | which profiles to run |
| `--no-build` | | reuse the existing build tree |
| `--cuda-arch <a>` | `native` | CUDA architecture |

## What gets measured

An instance is defined exactly as in the request: an initial direction, a final
direction, and N points in X and Y, plus the curvature and the DP settings
(`discr`, `nref`). All of it is written to `problems_<profile>.csv` **once** and
every solver reads that same file, so no two solvers can ever see different
inputs.

Each solver writes `results_<profile>_<solver>.csv` with, per instance:

- `length` — the length the solver's own DP bookkeeping reports;
- `length_check` — the length recomputed in double precision from the returned
  angles. This is the honest number: it is what the returned path actually
  measures. A gap between the two columns means the solver's internal
  accounting disagrees with its own answer;
- `time_ms` — wall-clock time of the solve call, including host-side work and,
  for the GPU solvers, all transfers;
- `time_solver_ms` — for the GPU solver, device time only (kernels plus
  transfers, excluding host-side angle sampling);
- `angles` — the N angles of the solution.

Instances not reached before `--max-seconds` are recorded as `skipped` rather
than dropped, so a slow solver's coverage is explicit.

### Profiles

**`budget`** sizes each configuration from a total curve budget. The DP
evaluates roughly `rounds * (2*discr + (N-3)*discr^2)` point-to-point curves per
instance — the first and last angles are fixed, so those two stages cost only
`discr` each — and the generator inverts that to pick `discr` per configuration.
With the default `--budget 1000000` the whole sweep costs about 1e6 curves.

**`stress`** keeps realistic discretisations (90 to 360) over the same N ladder
(3, 5, 10, 25, 50, 100, 200). It costs about 2.9e7 curves, still well inside the
60 s cap, and it is the profile where the GPU has enough work per launch to be
worth using. The `budget` profile's high-N entries end up at `discr = 12`, which
is far too little work to fill a GPU — that is a property of a 1e6 budget spread
over seven configurations, not of the implementations.

## Results on this machine

RTX 4070 Ti (sm_89, 60 SMs), CUDA 12.6, single-threaded CPU baseline (the
`#pragma omp parallel for` in `DP::solveDPInner` does nothing unless the library
is built with `-DMPDP_OPENMP=ON`).

`stress` profile, seven instances, 2.9e7 curves:

| solver | total time | speed-up | worst rel. error vs CPU |
|---|---|---|---|
| `cpu` (srcCC) | 4.965 s | 1.0x | — |
| `gpu_fp64` | 0.111 s | 44.2x | 4.0e-16 |
| `gpu_fp32` | 0.012 s | 393.8x | 4.7e-07 |

The now-removed legacy GPU solver (`DP::solveDP` in the old `dp.cuh`) measured
4.460 s / 1.1x for its `solveDPMatrixAllocator` variant and 0.236 s / 21.1x for
`solveDPAllIn1`, both after the compile and `cudaMemcpy` fixes needed to make it
run at all. It was deleted once the comparison was done.

Per-instance, the fp64 solver runs 40x to 60x faster than the CPU; the
largest case (N=200, discr=90) goes from 1382 ms to 39 ms.

`fp32` reaches 394x but its DP accumulates in single precision, so its reported
`length` drifts from the true path length by up to 1.8e-2 on the largest
instance. The path it picks is still within 4.7e-7 relative of the CPU's, so the
useful pattern is to search in fp32 and take `length_check` (the double-precision
recomputation over the returned angles) as the answer.

## Using the GPU solver

```cpp
#include <dp.cuh>

std::vector<Configuration2> points = /* N points; angle 0 where free */;
std::vector<bool> fixedAngles(points.size(), false);
fixedAngles.front() = fixedAngles.back() = true;   // initial/final direction
std::vector<real_type> params = {kmax};

mpdp::gpu::Options o;
o.discr     = 180;                              // angle samples per point
o.nref      = 4;                                // refinement rounds
o.precision = mpdp::gpu::Precision::FP64;       // or FP32
o.curve     = mpdp::gpu::CurveKind::DUBINS;

mpdp::gpu::Result r = mpdp::gpu::solveDP(points, fixedAngles, params, o);
// r.length, r.angles, r.device_ms, r.curves
```

Call `cudaFree(0)` once at startup so CUDA context creation (~50-200 ms) is not
charged to your first solve. `solveDP` throws `std::runtime_error` on a CUDA
failure or an infeasible problem; it allocates and frees its device buffers per
call, so if you solve many small instances in a loop that allocation shows up in
the profile.

### What actually drives the speed-up

Measured on the RTX 4070 Ti with the two variables separated:

| N (discr fixed at 180) | cpu | fp64 | fp32 |
|---|---|---|---|
| 5 | 60.6 ms | 42.6x | 171x |
| 25 | 617 ms | 58.0x | 634x |
| 200 | 5359 ms | 57.9x | 813x |

| discr (N fixed at 25) | cpu | fp64 | fp32 |
|---|---|---|---|
| 32 | 24.3 ms | 7.6x | 29x |
| 90 | 161 ms | 35.1x | 144x |
| 180 | 614 ms | 59.3x | 628x |
| 720 | 9600 ms | 69.5x | 2077x |

The fp64 speed-up is essentially **flat in N** (43x to 58x) and rises sharply
with **discr** (8x to 70x). A stage of the backward sweep costs `discr^2`
independent curve evaluations and is one kernel launch, so `discr` sets how much
parallelism each launch has, while N only adds more sequential launches. Use
`discr >= 90`; below about 64 the grid is too small to fill the card and the
per-launch overhead dominates.

N = 3 is the degenerate case: with both endpoints fixed there is no `discr^2`
stage at all, only two `discr`-sized ones, and the GPU is barely ahead (1.4x).

fp32 keeps improving with N as well, because at 1-6 ms per solve the fixed
per-launch cost is a larger share and N amortises it.

### Choosing the precision

Use **FP64** when the length you report has to be the CPU's, or when the
instance is coarse (`discr < 64`, `nref <= 1`) — there, near-ties decide the
answer and single precision flips them more often.

Use **FP32** for search, planning loops, or dataset generation, and take the
length from a double-precision recomputation over the returned angles rather
than from `Result.length`:

```cpp
double len = 0.0;
for (size_t i = 0; i + 1 < points.size(); ++i)
    len += mpdp::gpu::Dubins<double>::length(
        points[i].x(), points[i].y(), r.angles[i],
        points[i+1].x(), points[i+1].y(), r.angles[i+1], params.data());
```

That costs O(N) and removes the only real fp32 drawback: the DP accumulates the
running cost in float, so `Result.length` drifts (up to 1.8e-2 on N=200) even
though the *path* it selects is within ~5e-7 relative of the CPU's.

## Reeds-Shepp

`srcCU/include/rs.cuh` ports `RS::reeds_shepp` from `srcCC/rs.cc`: twelve base
words, each tried with four sign variants of the standardised problem, giving
the classic 48 maneuvers. The dispatch order and maneuver numbering match the
CPU exactly.

Select it with `mpdp::gpu::CurveKind::REEDS_SHEPP`, or from the demo:

```bash
./build/MPDPCU_exec --curve rs --discr 180
```

`MPMDBenchmarkRSCheck` validates the port point to point against the CPU `RS`:

| | result over 200000 random configurations |
|---|---|
| fp64 length vs CPU | **bit-identical** (worst relative difference 0.0) |
| fp64 maneuver number vs CPU | identical on every sample |
| maneuvers exercised | 46 of 48 |
| fp32 length vs CPU | worst relative difference 2.1e-04 |
| fp32 maneuver number vs CPU | differs on 0.09% of samples |

Two things worth knowing:

- **Only 46 maneuvers are reachable, on the CPU as well.** Maneuvers 2 and 4 are
  switched off in `srcCC/rs.cc` by `if (var < length && false)` — the `&& false`
  looks like a debugging leftover. The GPU keeps them disabled (the code is
  present but commented, next to an explanation) so the two agree; re-enabling
  them is a one-line change on each side. Worth a look, since it may mean the
  CPU is missing some optimal C|C|C paths.
- **Reeds-Shepp in fp32 is far less trustworthy than Dubins in fp32** (2.1e-04
  against ~5e-07). The base words branch on tight tolerances — `EPS3 = 1e-14`,
  and a `fabs(va) < 0.001` guard in `cc_c` — and single precision flips those
  decisions. Use fp64 for Reeds-Shepp unless you have measured that it does not
  matter for your case.

On the GPU, Reeds-Shepp costs roughly 6x a Dubins solve on the same problem
(48 maneuvers instead of 6): the Kaya 4 example at `discr = 180` takes 6.0 ms
with Dubins and 38.0 ms with Reeds-Shepp, and returns a shorter path (6.596
against 7.468) because it may reverse.

### No multi-point CPU baseline for Reeds-Shepp

The comparison above is point to point, because `DP::solveDP<RS>` does not
compile on the CPU: `solveDPInner` builds curves as
`CurveT(x0, y0, th0, x1, y1, th1, params)` and `RS` only offers the
`(Configuration2, Configuration2, params)` constructor. Adding the seven
argument constructor to `srcCC/include/rs.hh` would be enough to get a CPU
baseline for the multi-point case; it has not been done here.

The correctness argument in the meantime: the DP machinery is curve agnostic and
is validated against the CPU with Dubins, and the Reeds-Shepp curve is validated
bit-exactly point to point, so the composition of the two is sound.

## Angle sets, and the CPU bug this harness exposed

The GPU reproduces the CPU's sampling algorithm: the uniform samples, the
per-pair "guess" angles (the heading that lines two points up, and the tangents
to the two circles of radius 1/Kmax through both), and the same shrinking
refinement window. Rows therefore have different lengths, so the DP matrix is
stored padded to the longest row with an explicit per-row length; kernels never
read past a row's length, so the padding costs memory but no work. Ties in the
reduction keep the smaller index, which is what the CPU's ascending scan does.

Building the GPU against the CPU this way turned up a bug in the CPU, **since
fixed**. `srcCC/include/dp.hh` used to declare

```cpp
static K_T Kmax = DUBINS_DEFAULT_KMAX;   // 0.01
```

at namespace scope in a header. `static` gives it internal linkage, so every
translation unit got its own copy. `DP::solveDP` is a template defined in the
header, so `Kmax = params[0]` wrote the copy belonging to whichever TU
instantiated it — while `DP::guessInitialAngles`, compiled into `dp.cc`, read
`dp.cc`'s own copy, which nothing ever assigned and which stayed at `0.01`. The
guess angles were therefore always built for a turning radius of 100, whatever
curvature was requested.

The effect was not cosmetic. With `kmax = 3` the circles of radius 1/3 through
two points 8 apart do not exist, so the correct algorithm contributes only the
straight-line angle and rows hold `discr + 2` entries. With the stale radius of
100 the circles always existed and rows held `discr + 10` — eight extra
near-straight-line candidates per row, which made the CPU both slower and, at
coarse settings, accidentally slightly better than the algorithm as designed.

The diagnosis was confirmed against an independent third implementation: a plain
Python DP reproduced the GPU's answer bit-for-bit using `kmax`, and the old
CPU's answer bit-for-bit using `0.01`.

`Kmax` is now the `DP::k_max_` member. After the fix, over 12 instances spanning
N in 5..50 and `discr` in 32..180:

| nref | CPU before fix vs GPU fp64 | CPU after fix vs GPU fp64 |
|---|---|---|
| 0 | 9/12 exact, worst 4.8e-07 | **12/12 bit-identical** |
| 1 | 9/12 exact, worst 2.0e-08 | **12/12 bit-identical** |
| 2 | 10/12 exact, worst 9.5e-11 | **12/12 bit-identical** |
| 4 | 7/12 exact, worst 5.3e-16 | 7/12 exact, worst 4.3e-16 |

With no refinement the two implementations are now identical to the last bit.
The residual at `nref = 4` is ordinary floating-point noise: CUDA's `sin`,
`cos`, `atan2` and `acos` are not bit-identical to glibc's, and five rounds of a
discrete argmin occasionally latch onto a different member of a near-tie. (It is
not FMA contraction — building with `-fmad=false` changes nothing.)

Impact of the fix on the CPU's own results: at `discr >= 90` every instance in
both profiles returns an **identical** length, and the CPU got 1.03x to 1.18x
faster. Only at coarse discretisations does it differ at all — one instance
(N=25, `discr=36`) came out 5.15e-09 relative longer, having lost the accidental
extra candidates — and there the CPU is up to 2.4x faster.

Checked independently against the published Kaya/Omega optima (`examples/MPMD`),
20 configurations spanning `discr` 16/90 and `nref` 4/8:

- 8 identical, and every configuration 1.0x to 1.75x faster;
- **Kaya 4 improved substantially and consistently**: the error against the
  published optimum drops from 3.86e-04 to 2.02e-06 at `discr = 90` (both
  `nref` values), and from 3.95e-03 to 4.93e-04 at `discr = 16`;
- the apparent regressions on Kaya 2 and 3 are at 1e-15, i.e. machine epsilon;
- Omega moves from 7.53e-04 to 1.24e-03 at `discr = 90`, but that is refinement
  jitter rather than degradation: its error is non-monotonic in `discr` for
  *both* versions (post-fix 1.24e-03, 2.56e-04, 1.27e-03, 1.95e-04 at
  `discr` = 90, 100, 120, 180; pre-fix 7.53e-04, 2.56e-04, 1.28e-03, 1.70e-04
  over the same ladder). The residual is set by which local optimum the greedy
  refinement lands in, not by the sampling fix.

A second bug was fixed alongside it: `solveDP` called the first round as
`solveDPInner (compPoints, params)` without the template argument, so the coarse
round always ran on `Dubins` regardless of `CurveT` while the refinements
honoured it. Harmless while only Dubins was in use (nothing in the tree
instantiates `solveDP` with another curve today), fatal for any other curve.

## Files

| file | role |
|---|---|
| `run_benchmark.sh` | builds and runs everything (start here) |
| `gen_problems.cc` | writes the shared problem set |
| `bench_cpu.cc` | driver for `srcCC`'s `DP::solveDP` |
| `bench_gpu.cu` | driver for `mpdp::gpu::solveDP` (`--precision fp32` / `fp64`) |
| `rs_check.cu` | point-to-point validation of the GPU Reeds-Shepp port |
| `bench_io.hh` | problem/result structs and CSV I/O, no MPDP dependency |
| `bench_driver.hh` | shared main-loop, time cap, result writing |
| `summarize.py` | merges the result files into `summary_<profile>.md` |
| `run_benchmark.sh` | builds and runs everything |

The drivers are separate executables on purpose. Originally they had to be:
`srcCC` and `srcCU` both defined `Dubins`, `Curve` and `Configuration2` in the
global namespace, so linking both libraries into one binary was a
duplicate-symbol trap. That is fixed - the two trees now live in `mpdp::cpu` and
`mpdp::gpu` and do link together - but separate processes are still the right
shape here: CUDA context creation stays out of the CPU timings, and a crash in
one solver leaves the other's results intact.

Note that the repository's `.gitignore` excludes `*.csv`, so the generated data
files under `results/` are not tracked; the `summary_*.md` reports are.

## Namespaces

Everything is under `mpdp`:

| namespace | holds |
|---|---|
| `mpdp` | the shared typedefs (`real_type`, `Angle`, `LEN_T`, `K_T`, `ANGLE`), `TimePerf`, `AsyPlot`, the IO helpers |
| `mpdp::cpu` | all of `srcCC`: `Configuration2`, `Curve`, `Dubins`, `RS`, `DP`, the math utilities |
| `mpdp::gpu` | all of `srcCU`: `Configuration2`, `Dubins<T>`, `ReedsShepp<T>`, `Scalar<T>`, `solveDP` |

Because `mpdp::gpu` is nested inside `mpdp`, the GPU code still refers to
`real_type`, `Angle` and friends unqualified.

The examples, tests and executables open the namespaces they need with
`using namespace mpdp; using namespace mpdp::cpu;` (or `mpdp::gpu`) at the top,
so the demo code itself is unchanged.

This is an API break for anything consuming the exported `MPDPCC` target:
`Dubins` becomes `mpdp::cpu::Dubins`, and so on.

One duplicate remains: `srcCU/asyplot.cu` is a verbatim copy of
`srcCC/asyplot.cc`, so both libraries define `mpdp::AsyPlot`. It only bites a
program that links both libraries *and* uses `AsyPlot` - nothing does today.
Deleting `srcCU/asyplot.cu` and letting the CUDA side borrow the drawing code
from `MPDPCC` would close it.

## srcCU layout


The GPU tree now mirrors `srcCC` file for file:

| `srcCU` | mirrors | holds |
|---|---|---|
| `dp.cu` / `include/dp.cuh` | `dp.cc` / `dp.hh` | the DP kernels and `mpdp::gpu::solveDP` |
| `include/dubins.cuh` | `dubins.hh` | `mpdp::gpu::Dubins<T>` |
| `include/rs.cuh` | `rs.hh` | `mpdp::gpu::ReedsShepp<T>` |
| `include/curve.cuh` | `curve.hh` | `CurveKind` and the interface a curve family implements |
| `include/math_utils.cuh` | `math_utils.hh` | `Scalar<T>` and `mod2pi` |
| `include/configuration.cuh` | `configuration.hh` | `Configuration2` |

Everything a kernel calls lives in a header and is inlined, so the library no
longer needs relocatable device code (`CUDA_SEPARABLE_COMPILATION` now defaults
to `OFF`).
