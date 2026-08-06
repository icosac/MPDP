#!/usr/bin/env bash
#
# Builds and runs the CPU vs GPU comparison.
#
# Every file this produces lands in examples/MPMDBenchmark/results/. The only
# thing written elsewhere is the CMake build tree, which defaults to
# <repo>/build-mpmd-benchmark and can be pointed anywhere with BUILD_DIR.
#
# Usage:
#   ./examples/MPMDBenchmark/run_benchmark.sh [options]
#
# Options:
#   --budget <n>        Curve budget for the `budget` profile (default 1000000)
#   --max-seconds <s>   Wall-clock cap per solver per profile (default 60)
#   --instances <n>     Random instances per configuration (default 1)
#   --reps <n>          Timing repetitions per instance, fastest wins (default 3)
#   --profiles <list>   Comma separated: budget,stress (default both)
#   --no-build          Reuse the existing build tree
#   --cuda-arch <a>     CUDA architecture, e.g. 89 or native (default native)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
RESULTS_DIR="${SCRIPT_DIR}/results"
BUILD_DIR="${BUILD_DIR:-${REPO_DIR}/build-mpmd-benchmark}"

BUDGET=1000000
MAX_SECONDS=60
INSTANCES=1
REPS=3
PROFILES="budget,stress"
DO_BUILD=1
CUDA_ARCH="${CUDA_ARCH:-native}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --budget) BUDGET="$2"; shift 2 ;;
    --max-seconds) MAX_SECONDS="$2"; shift 2 ;;
    --instances) INSTANCES="$2"; shift 2 ;;
    --reps) REPS="$2"; shift 2 ;;
    --profiles) PROFILES="$2"; shift 2 ;;
    --no-build) DO_BUILD=0; shift ;;
    --cuda-arch) CUDA_ARCH="$2"; shift 2 ;;
    -h|--help) sed -n '2,25p' "$0"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

mkdir -p "${RESULTS_DIR}"

# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------
if [[ "${DO_BUILD}" -eq 1 ]]; then
  echo "=== Configuring in ${BUILD_DIR} ==="

  # /usr/bin/nvcc is often a symlink into /usr/local/cuda-*/bin, and nvcc
  # resolves its bundled headers relative to argv[0]; going through the symlink
  # can make it fail to find cuda_runtime.h. Prefer the real path.
  NVCC_ARG=()
  if command -v nvcc >/dev/null 2>&1; then
    NVCC_REAL="$(readlink -f "$(command -v nvcc)")"
    NVCC_ARG=(-DCMAKE_CUDA_COMPILER="${NVCC_REAL}")
  fi

  cmake -B "${BUILD_DIR}" -S "${REPO_DIR}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS="-O3 -std=c++20" \
    -DCUDA_ON=ON \
    -DCUDA_ARCH="${CUDA_ARCH}" \
    -DCOMPILE_EXAMPLES=ON \
    -DCOMPILE_CXX_EXEC=OFF \
    "${NVCC_ARG[@]}"

  echo "=== Building ==="
  cmake --build "${BUILD_DIR}" --target MPMDBenchmark -j "$(nproc)"
fi

GEN="${BUILD_DIR}/MPMDBenchmarkGen"
CPU="${BUILD_DIR}/MPMDBenchmarkCPU"
GPU="${BUILD_DIR}/MPMDBenchmarkGPU"

for exe in "${GEN}" "${CPU}"; do
  [[ -x "${exe}" ]] || { echo "Missing ${exe}; build first." >&2; exit 1; }
done

HAVE_GPU=1
if [[ ! -x "${GPU}" ]]; then
  echo "WARNING: GPU driver not built, running the CPU baseline only." >&2
  HAVE_GPU=0
fi

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
IFS=',' read -ra PROFILE_LIST <<< "${PROFILES}"
for profile in "${PROFILE_LIST[@]}"; do
  echo
  echo "############################################################"
  echo "### Profile: ${profile}"
  echo "############################################################"

  PROBLEMS="${RESULTS_DIR}/problems_${profile}.csv"
  "${GEN}" --out "${PROBLEMS}" --profile "${profile}" \
           --budget "${BUDGET}" --instances "${INSTANCES}"

  common=(--problems "${PROBLEMS}" --max-seconds "${MAX_SECONDS}" --reps "${REPS}")

  echo "--- CPU (srcCC) ---"
  "${CPU}" "${common[@]}" --out "${RESULTS_DIR}/results_${profile}_cpu.csv"

  if [[ "${HAVE_GPU}" -eq 1 ]]; then
    echo "--- GPU, fp64 ---"
    "${GPU}" "${common[@]}" --precision fp64 \
        --out "${RESULTS_DIR}/results_${profile}_gpu_fp64.csv"

    echo "--- GPU, fp32 ---"
    "${GPU}" "${common[@]}" --precision fp32 \
        --out "${RESULTS_DIR}/results_${profile}_gpu_fp32.csv"
  fi

  echo "--- Summary ---"
  python3 "${SCRIPT_DIR}/summarize.py" \
      --problems "${PROBLEMS}" --results "${RESULTS_DIR}" --profile "${profile}"
done

echo
echo "All results are under ${RESULTS_DIR}"
