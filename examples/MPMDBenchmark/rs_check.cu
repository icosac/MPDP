/**
 * @file rs_check.cu
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License 3.0.
 * @brief Point-to-point validation of the GPU Reeds-Shepp solver.
 *
 * Compares `mpdp::gpu::ReedsShepp<T>` against `srcCC`'s `RS` over random
 * configurations, in both precisions, and reports which of the 48 maneuvers
 * were exercised.
 *
 * This is the one place where the C++ and the CUDA trees are linked into the
 * same binary, and it only works because `rs.cuh` includes nothing but
 * `curve.cuh` and `math_utils.cuh` - in particular no `Configuration2`, which
 * is the type that clashes between the two trees. Compiled by nvcc but with no
 * device code: `ReedsShepp` is `__host__ __device__`, so the same source that
 * runs in the kernels runs here on the CPU.
 *
 * Usage: MPMDBenchmarkRSCheck [samples]   (default 200000)
 */
#include <cmath>
#include <cstdio>
#include <map>
#include <random>
#include <vector>

#include <configuration.hh>
#include <rs.hh>          // CPU
#include <rs.cuh>         // GPU (header-only, no Configuration2)

using namespace mpdp;
using namespace mpdp::cpu;

int main (int argc, char** argv)
{
  const int N = (argc > 1) ? std::atoi (argv[1]) : 200000;
  std::mt19937 gen (20240806);
  std::uniform_real_distribution<double> pos (-6.0, 6.0);
  std::uniform_real_distribution<double> ang (-M_PI, M_PI);
  std::uniform_real_distribution<double> kd (0.3, 5.0);

  int worst_man_mismatch = 0, len_mismatch = 0, both_inf = 0;
  double worst_rel = 0.0, worst_abs = 0.0;
  std::map<int,int> man_hist;
  double worst_case[7] = {0};
  double worst_rel_f = 0.0; int man_mismatch_f = 0;

  for (int i = 0; i < N; ++i)
  {
    double x0 = pos (gen), y0 = pos (gen), t0 = ang (gen);
    double x1 = pos (gen), y1 = pos (gen), t1 = ang (gen);
    double k  = kd (gen);

    RS cpu (Configuration2 (x0, y0, t0), Configuration2 (x1, y1, t1), {k});
    double lc = cpu.reeds_shepp (-1);

    const double params[1] = {k};
    auto g = mpdp::gpu::ReedsShepp<double>::solve (x0, y0, t0, x1, y1, t1, params);

    man_hist[g.man]++;

    const float fp[1] = {(float)k};
    auto gf = mpdp::gpu::ReedsShepp<float>::solve ((float)x0,(float)y0,(float)t0,
                                                   (float)x1,(float)y1,(float)t1, fp);
    if (std::isfinite (lc) && lc < 1e99 && std::isfinite (gf.len) && gf.len < 1e29)
    {
      double rf = std::fabs ((double)gf.len - lc) / (std::fabs (lc) > 1e-12 ? std::fabs (lc) : 1.0);
      if (rf > worst_rel_f) worst_rel_f = rf;
      if (gf.man != cpu.getNman()) ++man_mismatch_f;
    }

    const bool ci = !std::isfinite (lc) || lc > 1e99;
    const bool gi = !std::isfinite (g.len) || g.len > 1e99;
    if (ci && gi) { ++both_inf; continue; }

    if (cpu.getNman() != g.man) { ++worst_man_mismatch; }

    const double a = std::fabs (lc - g.len);
    const double r = a / (std::fabs (lc) > 1e-12 ? std::fabs (lc) : 1.0);
    if (a > 1e-9) { ++len_mismatch; }
    if (r > worst_rel)
    {
      worst_rel = r; worst_abs = a;
      worst_case[0]=x0; worst_case[1]=y0; worst_case[2]=t0;
      worst_case[3]=x1; worst_case[4]=y1; worst_case[5]=t1; worst_case[6]=k;
    }
  }

  printf ("samples              : %d\n", N);
  printf ("both infinite        : %d\n", both_inf);
  printf ("maneuver mismatches  : %d\n", worst_man_mismatch);
  printf ("length mismatches>1e-9: %d\n", len_mismatch);
  printf ("worst relative diff  : %.3e  (abs %.3e)\n", worst_rel, worst_abs);
  if (worst_rel > 0)
    printf ("  at ci=(%.6f,%.6f,%.6f) cf=(%.6f,%.6f,%.6f) k=%.6f\n",
            worst_case[0],worst_case[1],worst_case[2],
            worst_case[3],worst_case[4],worst_case[5],worst_case[6]);
  printf ("--- fp32 ---\n");
  printf ("worst relative diff  : %.3e\n", worst_rel_f);
  printf ("maneuver mismatches  : %d (%.2f%%)\n", man_mismatch_f, 100.0*man_mismatch_f/N);
  printf ("distinct maneuvers exercised: %zu -> ", man_hist.size());
  for (auto& kv : man_hist) printf ("%d ", kv.first);
  printf ("\n");
  return (worst_man_mismatch == 0 && len_mismatch == 0) ? 0 : 1;
}
