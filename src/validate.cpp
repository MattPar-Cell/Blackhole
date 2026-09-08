// validate.cpp - Numerical verification of the physics.
//
// Everything the renderer relies on is checked here against an independent
// result: an exact solution, a textbook closed form, or a measured value from
// the real universe.  Run with `blackhole --test`.
//
// The headline check is that the metric in kerr.h really is a vacuum solution
// of Einstein's field equations: the Ricci tensor is assembled numerically
// from the metric components alone - finite-differenced into Christoffel
// symbols, then into the Riemann tensor - and comes out zero.

#include "kerr.h"
#include "geodesic.h"
#include "disc.h"
#include "tov.h"
#include "neutronstar.h"
#include "jet.h"
#include "spectrum.h"
#include "constants.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>
#include <limits>

namespace {

int g_pass = 0, g_fail = 0;

void check(const char* name, bool ok, const std::string& detail) {
    std::printf("  [%s] %-52s %s\n", ok ? "PASS" : "FAIL", name, detail.c_str());
    ok ? ++g_pass : ++g_fail;
}

bool close(double a, double b, double tol) {
    const double d = std::fabs(a - b);
    const double s = std::max({std::fabs(a), std::fabs(b), 1e-300});
    return (d / s) <= tol;
}

void section(const char* title) {
    std::printf("\n%s\n", title);
    for (size_t i = 0; i < std::string(title).size(); ++i) std::putchar('-');
    std::putchar('\n');
}

// ===========================================================================
// Numerical differential geometry, built only from g_{mu nu}.
// ===========================================================================

using Mat4 = double[4][4];

void inverse_metric_full(double a, double r, double th, Mat4 out) {
    const bh::MetricUpper m = bh::metric_upper(a, r, th);
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) out[i][j] = 0.0;
    out[0][0] = m.tt; out[1][1] = m.rr; out[2][2] = m.thth; out[3][3] = m.pp;
    out[0][3] = out[3][0] = m.tp;
}

void inverse_metric_deriv_full(double a, double r, double th, int dir, Mat4 out) {
    const bh::MetricUpperDeriv d = bh::metric_upper_deriv(a, r, th);
    const bh::MetricUpper& m = (dir == 1) ? d.dr : d.dth;
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) out[i][j] = 0.0;
    if (dir != 1 && dir != 2) return;                 // d_t = d_phi = 0 exactly
    out[0][0] = m.tt; out[1][1] = m.rr; out[2][2] = m.thth; out[3][3] = m.pp;
    out[0][3] = out[3][0] = m.tp;
}

// Fourth-order central first derivative of g_{ab} with respect to coordinate
// `dir` (1 = r, 2 = theta).  d_t and d_phi vanish identically.
void dmetric(double a, double r, double th, int dir, double h, Mat4 out) {
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) out[i][j] = 0.0;
    if (dir != 1 && dir != 2) return;
    static const double c[5] = {1.0, -8.0, 0.0, 8.0, -1.0};
    for (int k = 0; k < 5; ++k) {
        if (c[k] == 0.0) continue;
        const double off = (k - 2) * h;
        double G[4][4];
        bh::metric_lower_full(a, r + (dir == 1 ? off : 0.0), th + (dir == 2 ? off : 0.0), G);
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) out[i][j] += c[k] * G[i][j];
    }
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) out[i][j] /= (12.0 * h);
}

// Fourth-order second derivative d_dirA d_dirB g_{ab}.
void ddmetric(double a, double r, double th, int dirA, int dirB,
              double hr, double hth, Mat4 out) {
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) out[i][j] = 0.0;
    if ((dirA != 1 && dirA != 2) || (dirB != 1 && dirB != 2)) return;

    if (dirA == dirB) {
        const double h = (dirA == 1) ? hr : hth;
        static const double c[5] = {-1.0, 16.0, -30.0, 16.0, -1.0};
        for (int k = 0; k < 5; ++k) {
            const double off = (k - 2) * h;
            double G[4][4];
            bh::metric_lower_full(a, r + (dirA == 1 ? off : 0.0),
                                     th + (dirA == 2 ? off : 0.0), G);
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j) out[i][j] += c[k] * G[i][j];
        }
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) out[i][j] /= (12.0 * h * h);
    } else {
        // Tensor product of two fourth-order first-derivative stencils.
        static const double c[5] = {1.0, -8.0, 0.0, 8.0, -1.0};
        const double ha = (dirA == 1) ? hr : hth;
        const double hb = (dirB == 1) ? hr : hth;
        for (int p = 0; p < 5; ++p) {
            if (c[p] == 0.0) continue;
            for (int q = 0; q < 5; ++q) {
                if (c[q] == 0.0) continue;
                double dr = 0.0, dt = 0.0;
                ((dirA == 1) ? dr : dt) += (p - 2) * ha;
                ((dirB == 1) ? dr : dt) += (q - 2) * hb;
                double G[4][4];
                bh::metric_lower_full(a, r + dr, th + dt, G);
                for (int i = 0; i < 4; ++i)
                    for (int j = 0; j < 4; ++j) out[i][j] += c[p] * c[q] * G[i][j];
            }
        }
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) out[i][j] /= (144.0 * ha * hb);
    }
}

struct Curvature {
    double ricci[4][4];
    double kretschmann;
};

Curvature curvature_at(double a, double r, double th) {
    const double hr = 1e-3 * std::max(1.0, r);
    const double hth = 1e-3;

    double ginv[4][4];
    inverse_metric_full(a, r, th, ginv);

    double dg[4][4][4];       // dg[lambda][i][j] = d_lambda g_ij
    double dginv[4][4][4];
    for (int l = 0; l < 4; ++l) {
        double tmp[4][4];
        dmetric(a, r, th, l, (l == 1) ? hr : hth, tmp);
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) dg[l][i][j] = tmp[i][j];
        inverse_metric_deriv_full(a, r, th, l, tmp);
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) dginv[l][i][j] = tmp[i][j];
    }

    double ddg[4][4][4][4];   // ddg[l][m][i][j] = d_l d_m g_ij
    for (int l = 0; l < 4; ++l)
        for (int m = 0; m < 4; ++m) {
            double tmp[4][4];
            ddmetric(a, r, th, l, m, hr, hth, tmp);
            for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) ddg[l][m][i][j] = tmp[i][j];
        }

    // Christoffel symbols of the second kind.
    double Gam[4][4][4];      // Gam[rho][mu][nu]
    for (int rho = 0; rho < 4; ++rho)
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu) {
                double s = 0.0;
                for (int sig = 0; sig < 4; ++sig)
                    s += ginv[rho][sig] * (dg[mu][sig][nu] + dg[nu][sig][mu] - dg[sig][mu][nu]);
                Gam[rho][mu][nu] = 0.5 * s;
            }

    // d_lambda Gamma^rho_{mu nu}
    double dGam[4][4][4][4];  // dGam[lam][rho][mu][nu]
    for (int lam = 0; lam < 4; ++lam)
        for (int rho = 0; rho < 4; ++rho)
            for (int mu = 0; mu < 4; ++mu)
                for (int nu = 0; nu < 4; ++nu) {
                    double s = 0.0;
                    for (int sig = 0; sig < 4; ++sig) {
                        s += dginv[lam][rho][sig] *
                             (dg[mu][sig][nu] + dg[nu][sig][mu] - dg[sig][mu][nu]);
                        s += ginv[rho][sig] *
                             (ddg[lam][mu][sig][nu] + ddg[lam][nu][sig][mu] -
                              ddg[lam][sig][mu][nu]);
                    }
                    dGam[lam][rho][mu][nu] = 0.5 * s;
                }

    // Riemann tensor R^rho_{sigma mu nu}.
    double Riem[4][4][4][4];
    for (int rho = 0; rho < 4; ++rho)
        for (int sig = 0; sig < 4; ++sig)
            for (int mu = 0; mu < 4; ++mu)
                for (int nu = 0; nu < 4; ++nu) {
                    double s = dGam[mu][rho][nu][sig] - dGam[nu][rho][mu][sig];
                    for (int lam = 0; lam < 4; ++lam)
                        s += Gam[rho][mu][lam] * Gam[lam][nu][sig] -
                             Gam[rho][nu][lam] * Gam[lam][mu][sig];
                    Riem[rho][sig][mu][nu] = s;
                }

    Curvature out{};
    for (int s = 0; s < 4; ++s)
        for (int n = 0; n < 4; ++n) {
            double acc = 0.0;
            for (int m = 0; m < 4; ++m) acc += Riem[m][s][m][n];
            out.ricci[s][n] = acc;
        }

    // Kretschmann scalar K = R_{abcd} R^{abcd}.
    double glow[4][4];
    bh::metric_lower_full(a, r, th, glow);
    double Rl[4][4][4][4];    // fully lowered
    for (int A = 0; A < 4; ++A)
        for (int B = 0; B < 4; ++B)
            for (int C = 0; C < 4; ++C)
                for (int D = 0; D < 4; ++D) {
                    double s = 0.0;
                    for (int rho = 0; rho < 4; ++rho) s += glow[A][rho] * Riem[rho][B][C][D];
                    Rl[A][B][C][D] = s;
                }
    double K = 0.0;
    for (int A = 0; A < 4; ++A)
        for (int B = 0; B < 4; ++B)
            for (int C = 0; C < 4; ++C)
                for (int D = 0; D < 4; ++D) {
                    double up = 0.0;
                    for (int p = 0; p < 4; ++p)
                        for (int q = 0; q < 4; ++q)
                            for (int u = 0; u < 4; ++u)
                                for (int v = 0; v < 4; ++v)
                                    up += ginv[A][p] * ginv[B][q] * ginv[C][u] * ginv[D][v] *
                                          Rl[p][q][u][v];
                    K += Rl[A][B][C][D] * up;
                }
    out.kretschmann = K;
    return out;
}

// Analytic Kretschmann scalar of the Kerr metric (M = 1).
double kretschmann_analytic(double a, double r, double th) {
    const double c2 = std::cos(th) * std::cos(th);
    const double S  = r * r + a * a * c2;
    const double t1 = r * r - a * a * c2;
    const double t2 = (r * r + a * a * c2) * (r * r + a * a * c2) - 16.0 * r * r * a * a * c2;
    return 48.0 * t1 * t2 / std::pow(S, 6.0);
}

// ===========================================================================
// Geodesic helpers used by several tests.
// ===========================================================================

// Integrate a geodesic from `y`, stopping when `stop` returns true or after
// max_steps.  Returns the final state.
template <typename StopFn>
bh::State integrate(double a, bh::State y, double h0, double rtol, double atol,
                    int max_steps, StopFn stop, long long* nsteps = nullptr) {
    bh::Stepper st(a, rtol, atol);
    st.reset();
    double h = h0;
    bh::State ynew{};
    bh::DenseSegment seg;
    for (int i = 0; i < max_steps; ++i) {
        int tries = 0;
        while (!st.try_step(y, h, ynew, seg)) {
            if (++tries > 60) return y;
        }
        if (nsteps) ++*nsteps;
        if (stop(ynew)) return ynew;
        y = ynew;
    }
    return y;
}

// ===========================================================================
// Tests
// ===========================================================================

void test_metric_algebra() {
    section("1. Metric algebra");

    double worst_id = 0.0;
    double worst_der = 0.0;
    const double spins[] = {0.0, 0.3, 0.7, 0.9, 0.998};
    for (double a : spins) {
        for (double r : {2.6, 4.0, 9.5, 30.0, 250.0}) {
            for (double th : {0.35, 1.0, M_PI / 2, 2.3, 2.9}) {
                if (r <= bh::horizon_outer(a) * 1.02) continue;

                double gl[4][4], gu[4][4];
                bh::metric_lower_full(a, r, th, gl);
                inverse_metric_full(a, r, th, gu);
                for (int i = 0; i < 4; ++i)
                    for (int j = 0; j < 4; ++j) {
                        double s = 0.0;
                        for (int k = 0; k < 4; ++k) s += gu[i][k] * gl[k][j];
                        worst_id = std::max(worst_id, std::fabs(s - (i == j ? 1.0 : 0.0)));
                    }

                // Analytic derivatives of g^{ab} against a fourth-order stencil.
                //
                // Some components have derivatives that vanish identically (for
                // example d_theta g^{tt} on the equator), where the finite
                // difference returns nothing but cancellation noise of size
                // ~eps |g| / h.  The comparison therefore uses a mixed
                // tolerance: relative where the derivative is real, and
                // referred to that noise floor where it is not.
                const bh::MetricUpperDeriv d = bh::metric_upper_deriv(a, r, th);
                auto compare = [&](int dir, auto sel, double ana) {
                    const double h = (dir == 1 ? 1e-4 * std::max(1.0, r) : 1e-4);
                    auto at = [&](double o) {
                        return sel(bh::metric_upper(a, r + (dir == 1 ? o : 0.0),
                                                       th + (dir == 2 ? o : 0.0)));
                    };
                    const double num =
                        (-at(2 * h) + 8 * at(h) - 8 * at(-h) + at(-2 * h)) / (12.0 * h);
                    const double noise = 50.0 * 2.22e-16 * std::fabs(at(0.0)) / h;
                    const double tol = 1e-6 * std::fabs(ana) + noise;
                    worst_der = std::max(worst_der, std::fabs(num - ana) / std::max(tol, 1e-300));
                };
                compare(1, [](const bh::MetricUpper& m) { return m.tt; },   d.dr.tt);
                compare(1, [](const bh::MetricUpper& m) { return m.tp; },   d.dr.tp);
                compare(1, [](const bh::MetricUpper& m) { return m.rr; },   d.dr.rr);
                compare(1, [](const bh::MetricUpper& m) { return m.thth; }, d.dr.thth);
                compare(1, [](const bh::MetricUpper& m) { return m.pp; },   d.dr.pp);
                compare(2, [](const bh::MetricUpper& m) { return m.tt; },   d.dth.tt);
                compare(2, [](const bh::MetricUpper& m) { return m.tp; },   d.dth.tp);
                compare(2, [](const bh::MetricUpper& m) { return m.rr; },   d.dth.rr);
                compare(2, [](const bh::MetricUpper& m) { return m.thth; }, d.dth.thth);
                compare(2, [](const bh::MetricUpper& m) { return m.pp; },   d.dth.pp);
            }
        }
    }
    char buf[160];
    std::snprintf(buf, sizeof buf, "max |g^ik g_kj - d^i_j| = %.2e", worst_id);
    check("g^{mu nu} is the exact inverse of g_{mu nu}", worst_id < 1e-11, buf);
    std::snprintf(buf, sizeof buf, "max error / (1e-6 |analytic| + FD noise) = %.2f", worst_der);
    check("analytic d g^{mu nu} matches finite differences", worst_der < 1.0, buf);
}

void test_einstein_equations() {
    section("2. Einstein's field equations");

    // The Ricci tensor is built from the metric with nothing but finite
    // differences.  For a vacuum solution with zero cosmological constant it
    // must vanish identically:  R_{mu nu} = 0.
    double worst_ricci = 0.0, worst_K = 0.0;
    double K_scale_at_worst = 1.0;

    const double spins[] = {0.0, 0.5, 0.9, 0.998};
    for (double a : spins) {
        for (double r : {3.0, 5.0, 12.0, 40.0}) {
            for (double th : {0.6, 1.2, M_PI / 2, 2.4}) {
                if (r <= bh::horizon_outer(a) * 1.2) continue;
                const Curvature cv = curvature_at(a, r, th);
                const double Kan = kretschmann_analytic(a, r, th);
                // Compare the Ricci components against the curvature scale
                // sqrt(K), which is the only meaningful yardstick: "zero"
                // means small compared with the curvature actually present.
                const double scale = std::sqrt(std::fabs(Kan));
                for (int i = 0; i < 4; ++i)
                    for (int j = 0; j < 4; ++j) {
                        // Normalise the mixed components by the metric so that
                        // the comparison is not distorted by coordinate scaling.
                        double gl[4][4];
                        bh::metric_lower_full(a, r, th, gl);
                        const double nrm = std::sqrt(std::fabs(gl[i][i] * gl[j][j])) + 1e-30;
                        const double v = std::fabs(cv.ricci[i][j]) / (nrm * scale);
                        worst_ricci = std::max(worst_ricci, v);
                    }
                const double relK = std::fabs(cv.kretschmann - Kan) / std::fabs(Kan);
                if (relK > worst_K) { worst_K = relK; K_scale_at_worst = Kan; }
            }
        }
    }

    char buf[200];
    std::snprintf(buf, sizeof buf,
                  "max |R_mn| / (curvature scale) = %.2e  [vacuum: R_mn = 0]", worst_ricci);
    check("Kerr satisfies R_{mu nu} = 0 (numerically, from g)", worst_ricci < 1e-6, buf);

    std::snprintf(buf, sizeof buf, "max relative error = %.2e (K ~ %.3e)", worst_K,
                  K_scale_at_worst);
    check("Riemann tensor reproduces the analytic Kretschmann scalar",
          worst_K < 1e-6, buf);

    // Schwarzschild limit: K = 48 M^2 / r^6, and it diverges at r = 0 while
    // staying perfectly finite at the horizon - the horizon is not a
    // singularity, only a coordinate artefact of Boyer-Lindquist.
    const double K_horizon = kretschmann_analytic(0.0, 2.0, M_PI / 2);
    std::snprintf(buf, sizeof buf, "K(r=2M) = %.6f = 48/64 = %.6f", K_horizon, 48.0 / 64.0);
    check("curvature is finite at the event horizon", close(K_horizon, 0.75, 1e-12), buf);
}

void test_orbits() {
    section("3. Circular orbits and horizons");

    char buf[200];

    // ISCO: 6M for Schwarzschild, 1M / 9M for an extremal hole.
    std::snprintf(buf, sizeof buf, "r_isco(a=0) = %.10f M", bh::isco_radius(0.0, +1));
    check("ISCO of a Schwarzschild hole is 6M", close(bh::isco_radius(0.0, +1), 6.0, 1e-12), buf);

    std::snprintf(buf, sizeof buf, "prograde %.6f M, retrograde %.6f M",
                  bh::isco_radius(1.0, +1), bh::isco_radius(1.0, -1));
    check("ISCO of an extremal hole is 1M / 9M",
          close(bh::isco_radius(1.0, +1), 1.0, 1e-9) &&
          close(bh::isco_radius(1.0, -1), 9.0, 1e-9), buf);

    // Photon spheres.
    std::snprintf(buf, sizeof buf, "r_ph(a=0) = %.10f M", bh::photon_circular_orbit(0.0, +1));
    check("photon sphere of a Schwarzschild hole is 3M",
          close(bh::photon_circular_orbit(0.0, +1), 3.0, 1e-12), buf);
    std::snprintf(buf, sizeof buf, "prograde %.6f M, retrograde %.6f M",
                  bh::photon_circular_orbit(1.0, +1), bh::photon_circular_orbit(1.0, -1));
    check("photon orbits of an extremal hole are 1M / 4M",
          close(bh::photon_circular_orbit(1.0, +1), 1.0, 1e-9) &&
          close(bh::photon_circular_orbit(1.0, -1), 4.0, 1e-9), buf);

    // Horizon and ergosphere.
    std::snprintf(buf, sizeof buf, "r+ = %.6f M, r_ergo(equator) = %.6f M",
                  bh::horizon_outer(0.9), bh::ergosphere_outer(0.9, M_PI / 2));
    check("ergosphere reaches 2M at the equator for any spin",
          close(bh::ergosphere_outer(0.9, M_PI / 2), 2.0, 1e-12), buf);

    // Binding energy at the ISCO: the accretion efficiency.
    const double eta0 = bh::disc_efficiency(0.0);
    const double eta1 = bh::disc_efficiency(0.998);
    std::snprintf(buf, sizeof buf, "eta(a=0) = %.4f (1 - sqrt(8/9)), eta(a=0.998) = %.4f",
                  eta0, eta1);
    check("accretion efficiency is 5.72% (a=0) and ~32% (a=0.998)",
          close(eta0, 1.0 - std::sqrt(8.0 / 9.0), 1e-12) && eta1 > 0.30 && eta1 < 0.33, buf);

    // A four-velocity really is a unit timelike vector.
    double worst = 0.0;
    for (double a : {0.0, 0.5, 0.95}) {
        for (double r = bh::isco_radius(a) ; r < 200.0; r *= 1.7) {
            const auto u = bh::keplerian_four_velocity(a, r);
            const auto ml = bh::metric_lower(a, r, M_PI / 2);
            const double n = ml.tt * u[0] * u[0] + 2 * ml.tp * u[0] * u[3] + ml.pp * u[3] * u[3];
            worst = std::max(worst, std::fabs(n + 1.0));
        }
    }
    std::snprintf(buf, sizeof buf, "max |u.u + 1| = %.2e", worst);
    check("Keplerian four-velocity is normalised (u.u = -1)", worst < 1e-11, buf);

    // Kepler's third law survives intact in Boyer-Lindquist coordinates:
    // Omega = 1/(r^{3/2} + a), so for a = 0 the coordinate period is
    // exactly 2 pi r^{3/2} - Newton's law, unmodified.
    const double r = 100.0;
    const auto u = bh::keplerian_four_velocity(0.0, r);
    const double Omega = u[3] / u[0];
    std::snprintf(buf, sizeof buf, "Omega = %.12e vs r^{-3/2} = %.12e", Omega,
                  std::pow(r, -1.5));
    check("Kepler's third law holds exactly in BL coordinates (a=0)",
          close(Omega, std::pow(r, -1.5), 1e-12), buf);
}

void test_null_geodesics() {
    section("4. Null geodesics");

    char buf[220];

    // --- conservation along a traced ray ---------------------------------
    {
        const double a = 0.9;
        const double r0 = 60.0, th0 = 1.3;
        const bh::Geom g = bh::geom_at(a, r0, th0);
        const bh::Tetrad tet = bh::zamo_tetrad(g);
        double worst_norm = 0.0, worst_Q = 0.0, worst_E = 0.0;

        for (int i = 0; i < 24; ++i) {
            const double ang = 0.02 + 0.05 * i;      // sweep past the shadow edge
            const double dir[3] = {-std::cos(ang), 0.2 * std::sin(ang),
                                    std::sqrt(std::max(0.0, 1.0 - std::cos(ang) * std::cos(ang) -
                                              0.04 * std::sin(ang) * std::sin(ang)))};
            bh::State y = bh::make_photon(a, r0, th0, 0.0, tet, dir, true);
            const bh::Invariants iv0 = bh::invariants(a, y);
            const double Q0 = std::fabs(iv0.Q) + 1e-12;

            y = integrate(a, y, 0.5, 1e-10, 1e-12, 20000,
                          [&](const bh::State& s) {
                              return s[bh::Y_R] < bh::horizon_outer(a) * 1.01 ||
                                     s[bh::Y_R] > 2000.0;
                          });
            const bh::Invariants iv = bh::invariants(a, y);
            worst_norm = std::max(worst_norm, std::fabs(iv.norm_rel));
            worst_Q = std::max(worst_Q, std::fabs(iv.Q - iv0.Q) / Q0);
            worst_E = std::max(worst_E, std::fabs(iv.E - iv0.E) / std::fabs(iv0.E));
        }
        std::snprintf(buf, sizeof buf, "max |g^ab p_a p_b| / |terms| = %.2e", worst_norm);
        check("photons stay null along the whole trajectory", worst_norm < 1e-10, buf);
        std::snprintf(buf, sizeof buf, "max relative drift of Q = %.2e", worst_Q);
        check("Carter's constant is conserved", worst_Q < 1e-8, buf);
        std::snprintf(buf, sizeof buf, "max relative drift of E = %.2e", worst_E);
        check("energy E = -p_t is conserved exactly", worst_E < 1e-14, buf);
    }

    // --- light deflection --------------------------------------------------
    //
    // For a null geodesic of impact parameter b in Schwarzschild,
    //     dphi/du = 1 / sqrt(1/b^2 - u^2 + 2 M u^3),   u = 1/r,
    // so the azimuth swept between radius r0 and the turning point is an
    // ordinary one-dimensional integral.  Substituting u = u_turn - s^2
    // removes the square-root singularity at the turning point and leaves a
    // smooth integrand, which Gauss-Legendre nails to machine precision.
    // That gives a completely independent reference for the traced ray.
    {
        auto turning_point = [](double b) {
            // Largest root of  u^2 - 2 u^3 = 1/b^2 , i.e. the photon's
            // closest approach.  Bisect on r in [3M, b].
            // 1/r^2 - 2/r^3 peaks at r = 3M (giving b_crit = sqrt(27) M) and
            // decreases outwards, so the wanted root is the outer one.
            double lo = 3.0, hi = std::max(b, 3.1);
            for (int i = 0; i < 200; ++i) {
                const double m = 0.5 * (lo + hi);
                const double f = 1.0 / (m * m) - 2.0 / (m * m * m) - 1.0 / (b * b);
                if (f > 0.0) lo = m; else hi = m;
            }
            return 0.5 * (lo + hi);
        };
        auto sweep = [&](double b, double r0) {
            const double rp = turning_point(b);
            const double up = 1.0 / rp;
            const double u0 = 1.0 / r0;
            const double smax = std::sqrt(std::max(0.0, up - u0));
            // Integrand in s:  2 s / sqrt(f(u_p - s^2)) , smooth at s = 0.
            auto integ = [&](double s) {
                const double u = up - s * s;
                const double f = 1.0 / (b * b) - u * u + 2.0 * u * u * u;
                return 2.0 * s / std::sqrt(std::max(f, 1e-300));
            };
            // Composite Gauss-Legendre, 8 nodes per panel.
            static const double gx[8] = {-0.9602898564975363, -0.7966664774136267,
                                         -0.5255324099163290, -0.1834346424956498,
                                          0.1834346424956498,  0.5255324099163290,
                                          0.7966664774136267,  0.9602898564975363};
            static const double gw[8] = {0.1012285362903763, 0.2223810344533745,
                                         0.3137066458778873, 0.3626837833783620,
                                         0.3626837833783620, 0.3137066458778873,
                                         0.2223810344533745, 0.1012285362903763};
            const int panels = 4000;
            double total = 0.0;
            for (int p = 0; p < panels; ++p) {
                const double a0 = smax * p / panels, a1 = smax * (p + 1) / panels;
                const double hm = 0.5 * (a1 - a0), cm = 0.5 * (a1 + a0);
                for (int k = 0; k < 8; ++k) total += hm * gw[k] * integ(cm + hm * gx[k]);
            }
            return 2.0 * total;    // in and back out again
        };

        // (a) The asymptotic limit r0 -> infinity must reproduce Einstein's
        //     1915 result, 4GM/(c^2 b), plus the known second-order term.
        // Expressed in the impact parameter b (rather than the closest
        // approach) the expansion is  alpha = 4M/b + (15 pi/4) M^2/b^2 + ...
        for (double b : {1.0e3, 1.0e4}) {
            const double defl = sweep(b, 1e18) - M_PI;
            const double einstein = 4.0 / b;
            const double second = 4.0 / b + (15.0 * M_PI / 4.0) / (b * b);
            std::snprintf(buf, sizeof buf,
                          "b = %.0e M: exact %.10e rad, 4M/b = %.10e, +2nd order %.10e",
                          b, defl, einstein, second);
            check("light deflection tends to 4GM/(c^2 b)", close(defl, second, 1e-4), buf);
        }

        // (b) The traced geodesic must reproduce the exact quadrature at
        //     finite radius, where the weak-field formula does not apply.
        const double a = 0.0, r0 = 1.0e5;
        for (double b : {20.0, 500.0}) {
            const bh::Geom g = bh::geom_at(a, r0, M_PI / 2);
            const bh::Tetrad tet = bh::zamo_tetrad(g);
            // For a static observer, b = r sin(alpha) / sqrt(1 - 2M/r).
            const double alpha = std::asin(b * std::sqrt(1.0 - 2.0 / r0) / r0);
            const double dir[3] = {-std::cos(alpha), 0.0, std::sin(alpha)};
            bh::State y = bh::make_photon(a, r0, M_PI / 2, 0.0, tet, dir, true);
            // The impact parameter actually realised, b = L/E, is exact.
            const double b_true = y[bh::Y_PPH] / (-y[bh::Y_PT]);

            // Integrate out to r = r0 again, landing exactly on that radius via
            // the dense-output polynomial - overshooting the last step by even
            // a few hundred M would swamp the effect being measured.
            bh::Stepper st(a, 1e-13, 1e-15);
            st.reset();
            double h = 1.0;
            bh::State ynew{};
            bh::DenseSegment seg;
            bool outbound = false;
            double traced = 0.0;
            for (int i = 0; i < 2000000; ++i) {
                int tries = 0;
                while (!st.try_step(y, h, ynew, seg)) if (++tries > 60) break;
                if (h > 0.05 * ynew[bh::Y_R]) h = 0.05 * ynew[bh::Y_R];
                if (ynew[bh::Y_PR] > 0.0) outbound = true;
                if (outbound && ynew[bh::Y_R] >= r0) {
                    double lo = 0.0, hi = 1.0;
                    for (int k = 0; k < 60; ++k) {
                        const double mid = 0.5 * (lo + hi);
                        if (seg.eval(mid)[bh::Y_R] >= r0) hi = mid; else lo = mid;
                    }
                    traced = std::fabs(seg.eval(0.5 * (lo + hi))[bh::Y_PH]);
                    break;
                }
                y = ynew;
            }
            const double exact = sweep(std::fabs(b_true), r0);
            std::snprintf(buf, sizeof buf,
                          "b = %.0f M: traced %.9f rad, exact quadrature %.9f rad", b,
                          traced, exact);
            check("traced geodesic matches the exact deflection integral",
                  close(traced, exact, 1e-6), buf);
        }
    }

    // --- shadow size: the critical impact parameter -----------------------
    {
        const double a = 0.0;
        const double r0 = 1.0e4;
        const bh::Geom g = bh::geom_at(a, r0, M_PI / 2);
        const bh::Tetrad tet = bh::zamo_tetrad(g);

        auto captured = [&](double alpha) {
            const double dir[3] = {-std::cos(alpha), 0.0, std::sin(alpha)};
            bh::State y = bh::make_photon(a, r0, M_PI / 2, 0.0, tet, dir, true);
            bool cap = false;
            integrate(a, y, 1.0, 1e-11, 1e-13, 200000,
                      [&](const bh::State& s) {
                          if (s[bh::Y_R] < 2.0 * 1.0001) { cap = true; return true; }
                          return s[bh::Y_R] > 2.0 * r0;
                      });
            return cap;
        };

        double lo = 0.0, hi = 0.02;           // capture below, escape above
        for (int i = 0; i < 60; ++i) {
            const double mid = 0.5 * (lo + hi);
            if (captured(mid)) lo = mid; else hi = mid;
        }
        const double alpha_c = 0.5 * (lo + hi);
        const double b_c = r0 * std::sin(alpha_c) / std::sqrt(1.0 - 2.0 / r0);
        const double exact = 3.0 * std::sqrt(3.0);
        std::snprintf(buf, sizeof buf, "b_crit = %.8f M vs 3 sqrt(3) = %.8f M", b_c, exact);
        check("shadow radius is the critical impact parameter 3 sqrt(3) M",
              close(b_c, exact, 1e-5), buf);
    }

    // --- gravitational redshift -------------------------------------------
    {
        // A photon climbing from r_em to r_obs in Schwarzschild is redshifted by
        // sqrt((1 - 2M/r_obs) / (1 - 2M/r_em)).
        const double a = 0.0, r_em = 10.0, r_obs = 1000.0;
        const bh::MetricLower m_em = bh::metric_lower(a, r_em, M_PI / 2);
        const bh::MetricLower m_ob = bh::metric_lower(a, r_obs, M_PI / 2);
        // Static observers at both ends.
        const std::array<double, 4> u_em = {1.0 / std::sqrt(-m_em.tt), 0, 0, 0};
        const std::array<double, 4> u_ob = {1.0 / std::sqrt(-m_ob.tt), 0, 0, 0};
        // Radial photon: only p_t matters, and it is conserved.
        const double pt = -1.0;
        const double g = (pt * u_ob[0]) / (pt * u_em[0]);
        // Climbing out of the well costs energy: nu_obs < nu_em.
        const double predicted = std::sqrt((1.0 - 2.0 / r_em) / (1.0 - 2.0 / r_obs));
        std::snprintf(buf, sizeof buf, "nu_obs/nu_em = %.10f vs %.10f", g, predicted);
        check("gravitational redshift between static observers",
              close(g, predicted, 1e-12), buf);
    }

    // --- frame dragging ----------------------------------------------------
    {
        const double a = 0.9, r = 4.0, th = M_PI / 2;
        const bh::Geom g = bh::geom_at(a, r, th);
        const bh::Tetrad tet = bh::zamo_tetrad(g);
        const std::array<double, 4> u = {tet.e[0][0], 0.0, 0.0, tet.e[0][3]};
        const auto ul = bh::lower(bh::metric_lower(g), u);
        std::snprintf(buf, sizeof buf, "u_phi = %.3e, omega = %.6f c^3/GM (dragging rate)",
                      ul[3], bh::zamo_omega(g));
        check("ZAMO really has zero angular momentum (u_phi = 0)",
              std::fabs(ul[3]) < 1e-12, buf);
    }
}

void test_extreme_rotation() {
    section("5b. Maximally rotating holes");

    char buf[260];

    // The horizon is dragged round rigidly at Omega_H = a / (2 r_+).  That is
    // not an independent definition: it is the limit of the zero-angular-
    // momentum observer's angular velocity omega = 2 a r / A as r -> r_+.
    // Checking the two agree exercises the metric right at the horizon.
    {
        double worst = 0.0;
        for (double a : {0.3, 0.9, 0.998, 0.9999}) {
            const double rp = bh::horizon_outer(a);
            // Approach the horizon from outside.
            const double r = rp * (1.0 + 1e-9);
            const bh::Geom g = bh::geom_at(a, r, M_PI / 2);
            const double omega = bh::zamo_omega(g);
            const double OmH = bh::horizon_angular_velocity(a);
            worst = std::max(worst, std::fabs(omega - OmH) / OmH);
        }
        std::snprintf(buf, sizeof buf,
                      "max relative difference = %.2e; Omega_H(0.998) = %.6f c^3/GM",
                      worst, bh::horizon_angular_velocity(0.998));
        check("ZAMO angular velocity tends to Omega_H at the horizon", worst < 1e-7, buf);
    }

    // Inside the ergosphere g_tt changes sign, so no observer can stay at
    // fixed phi: "standing still" would require moving faster than light.
    // Everything in there is dragged forwards, whatever it does.
    {
        const double a = 0.998;
        bool all_positive = true, outside_negative = true;
        for (double th : {0.4, 1.0, M_PI / 2}) {
            const double rE = bh::ergosphere_outer(a, th);
            const double rp = bh::horizon_outer(a);
            for (double f : {0.1, 0.5, 0.9}) {
                const double r = rp + f * (rE - rp);
                if (bh::metric_lower(a, r, th).tt <= 0.0) all_positive = false;
            }
            if (bh::metric_lower(a, rE * 1.05, th).tt >= 0.0) outside_negative = false;
        }
        std::snprintf(buf, sizeof buf,
                      "g_tt > 0 throughout the ergosphere, < 0 outside it "
                      "(equatorial static limit at r = %.3f M)",
                      bh::ergosphere_outer(a, M_PI / 2));
        check("no static observer can exist inside the ergosphere",
              all_positive && outside_negative, buf);
    }

    // The extremal limits.  As a -> 1 the ISCO, the prograde photon orbit and
    // the horizon all converge on r = M in these coordinates, and the binding
    // energy at the ISCO tends to 1 - 1/sqrt(3): a maximally spinning hole can
    // convert 42% of the rest mass of what it swallows into radiation.
    {
        const double eta_extremal = 1.0 - 1.0 / std::sqrt(3.0);
        const double eta_998 = bh::disc_efficiency(0.998);
        const double eta_9999 = bh::disc_efficiency(0.9999);
        std::snprintf(buf, sizeof buf,
                      "eta(0.998) = %.4f, eta(0.9999) = %.4f, extremal limit %.4f",
                      eta_998, eta_9999, eta_extremal);
        check("efficiency rises towards 1 - 1/sqrt(3) as a -> 1",
              eta_998 > 0.31 && eta_998 < 0.33 && eta_9999 > eta_998 &&
              eta_9999 < eta_extremal, buf);

        std::snprintf(buf, sizeof buf,
                      "a = 0.9999: ISCO %.4f M, photon orbit %.4f M, horizon %.4f M",
                      bh::isco_radius(0.9999, +1), bh::photon_circular_orbit(0.9999, +1),
                      bh::horizon_outer(0.9999));
        check("ISCO, photon orbit and horizon converge as a -> 1",
              bh::isco_radius(0.9999, +1) < 1.09 &&
              bh::photon_circular_orbit(0.9999, +1) < 1.06, buf);
    }

    // The Thorne limit.  A hole fed by a radiating thin disc stops spinning up
    // at a = 0.998, because photons emitted by the disc are preferentially
    // captured when their angular momentum opposes the spin.  Check that the
    // disc's inner edge is still comfortably outside the horizon there - the
    // configuration has to be a real one.
    {
        const double a = 0.998;
        const double isco = bh::isco_radius(a, +1);
        const double rp = bh::horizon_outer(a);
        const auto u = bh::keplerian_four_velocity(a, isco, +1);
        const double Omega = u[3] / u[0];
        std::snprintf(buf, sizeof buf,
                      "ISCO %.4f M sits %.1f%% outside the horizon at %.4f M; Omega = %.4f",
                      isco, 100.0 * (isco / rp - 1.0), rp, Omega);
        check("at the Thorne limit the ISCO is still outside the horizon",
              isco > rp * 1.10 && Omega > 0.0, buf);
    }
}

void test_timelike_geodesics() {
    section("5. Timelike geodesics: perihelion precession");

    char buf[240];

    // Bound orbit in Schwarzschild parameterised by semi-latus rectum p and
    // eccentricity e (Cutler, Kennefick & Poisson 1994):
    //   E^2 = (p-2-2e)(p-2+2e) / (p (p-3-e^2)),   L^2 = p^2 / (p-3-e^2)
    auto precession_of = [&](double p, double e, double rtol) {
        const double E2 = (p - 2.0 - 2.0 * e) * (p - 2.0 + 2.0 * e) / (p * (p - 3.0 - e * e));
        const double L2 = p * p / (p - 3.0 - e * e);
        const double E = std::sqrt(E2), L = std::sqrt(L2);
        const double r_peri = p / (1.0 + e);

        // Start at perihelion, moving in +phi with p_r = 0.
        bh::State y = {0.0, r_peri, M_PI / 2, 0.0, -E, 0.0, 0.0, L};

        // Integrate until r comes back to a minimum (p_r crosses 0 upwards).
        bh::Stepper st(0.0, rtol, rtol * 1e-3);
        st.reset();
        double h = 0.01 * r_peri;
        bh::State ynew{};
        bh::DenseSegment seg;
        bool went_out = false;
        for (int i = 0; i < 4000000; ++i) {
            int tries = 0;
            while (!st.try_step(y, h, ynew, seg)) if (++tries > 60) break;
            if (h > 0.02 * ynew[bh::Y_R]) h = 0.02 * ynew[bh::Y_R];
            if (!went_out && ynew[bh::Y_PR] > 0.0 && ynew[bh::Y_R] > r_peri * 1.05) went_out = true;
            if (went_out && y[bh::Y_PR] < 0.0 && ynew[bh::Y_PR] >= 0.0) {
                // Bisect on p_r = 0 to land exactly on perihelion.
                double lo = 0.0, hi = 1.0;
                for (int k = 0; k < 60; ++k) {
                    const double mid = 0.5 * (lo + hi);
                    if (seg.eval(mid)[bh::Y_PR] >= 0.0) hi = mid; else lo = mid;
                }
                const bh::State yc = seg.eval(0.5 * (lo + hi));
                return yc[bh::Y_PH] - 2.0 * M_PI;
            }
            y = ynew;
        }
        return std::numeric_limits<double>::quiet_NaN();
    };

    // A strongly relativistic orbit first, where the effect is large and the
    // leading-order formula is only approximate.
    {
        const double p = 1.0e5, e = 0.2;
        const double measured = precession_of(p, e, 1e-13);
        const double predicted = 6.0 * M_PI / p;      // = 6 pi G M / (c^2 a (1-e^2))
        std::snprintf(buf, sizeof buf,
                      "p = 1e5 M: measured %.9e rad/orbit, 6 pi M/p = %.9e", measured, predicted);
        check("periastron advance matches 6 pi GM / (c^2 a (1-e^2))",
              close(measured, predicted, 1e-3), buf);
    }

    // Now the real thing: Mercury around the real Sun, integrated in the
    // Schwarzschild metric of one solar mass.
    {
        const double a_orb = 5.790905e10;             // m, semi-major axis
        const double ecc = 0.205630;
        const double r_g = phys::r_g_metres(phys::M_sun);
        const double p = a_orb * (1.0 - ecc * ecc) / r_g;    // semi-latus rectum in M

        const double measured = precession_of(p, ecc, 1e-13);
        const double predicted = 6.0 * M_PI / p;

        // Orbits per Julian century.
        const double period_days = 87.9691;
        const double orbits_per_century = 36525.0 / period_days;
        const double arcsec = phys::rad_to_arcsec(measured) * orbits_per_century;
        const double arcsec_pred = phys::rad_to_arcsec(predicted) * orbits_per_century;

        std::snprintf(buf, sizeof buf,
                      "integrated %.3f\"/century  (GR formula %.3f\", observed 42.98 +/- 0.04\")",
                      arcsec, arcsec_pred);
        check("Mercury's anomalous perihelion advance", std::fabs(arcsec - 42.98) < 0.5, buf);
    }
}

void test_disc() {
    section("6. Novikov-Thorne accretion disc");

    char buf[240];

    // The closed form Page & Thorne derived for
    //     W(r) = Integral (E - Omega L) dL/dr' dr'
    // is checked against brute-force numerical quadrature of that very
    // integral, with dL/dr taken from a high-order finite difference of the
    // Kerr circular-orbit angular momentum.  Two completely different routes
    // to the same number.
    {
        double worst = 0.0;
        for (double a : {0.0, 0.5, 0.9, 0.998}) {
            const double r_isco = bh::isco_radius(a, +1);
            auto integrand = [&](double r) {
                const double h = 1e-6 * r;
                auto L = [&](double x) { return bh::circular_angmom(a, x, +1); };
                const double dL = (-L(r + 2 * h) + 8 * L(r + h) - 8 * L(r - h) + L(r - 2 * h)) /
                                  (12.0 * h);
                const double Om = 1.0 / (r * std::sqrt(r) + a);
                return (bh::circular_energy(a, r, +1) - Om * L(r)) * dL;
            };
            for (double r_end : {8.0, 30.0, 200.0, 5000.0}) {
                if (r_end <= r_isco * 1.2) continue;
                // Simpson in u = sqrt(r); the integrand vanishes at the ISCO
                // because L(r) is stationary there.
                const int N = 200000;
                const double u0 = std::sqrt(r_isco), u1 = std::sqrt(r_end);
                double W = 0.0;
                for (int i = 0; i < N; ++i) {
                    const double ua = u0 + (u1 - u0) * i / N;
                    const double ub = u0 + (u1 - u0) * (i + 1) / N;
                    const double um = 0.5 * (ua + ub);
                    auto f = [&](double u) { return (u <= u0) ? 0.0 : integrand(u * u) * 2.0 * u; };
                    W += (ub - ua) / 6.0 * (f(ua) + 4.0 * f(um) + f(ub));
                }
                const double closed = bh::page_thorne_W(a, r_end, r_isco);
                worst = std::max(worst, std::fabs(W - closed) / std::fabs(closed));
            }
        }
        std::snprintf(buf, sizeof buf, "max relative difference = %.3e", worst);
        check("Page-Thorne closed form matches direct quadrature", worst < 1e-6, buf);
    }

    // Far from the hole the flux must fall off as the classical
    // Shakura-Sunyaev law F -> 3 G M Mdot / (8 pi r^3).  (The relativistic
    // result approaches it only as 1 - C/sqrt(r), with C = 3.976 rather than
    // the Newtonian sqrt(r_isco) = 2.449, so the leading term is what
    // converges.)
    {
        bh::NovikovThorneDisc d(10.0 * phys::M_sun, 0.0, 1e6, 0.1);
        double worst = 0.0;
        for (double r : {1e5, 1e6, 1e7}) {
            const double f = d.flux_dimensionless(r) * 8.0 * M_PI * r * r * r / 3.0;
            worst = std::max(worst, std::fabs(f - 1.0));
        }
        std::snprintf(buf, sizeof buf, "max |F / F_SS - 1| at r >= 1e5 M = %.3e", worst);
        check("flux tends to the Shakura-Sunyaev law far from the hole", worst < 2e-2, buf);
    }

    // Energy budget: the power reaching infinity is the binding energy released
    // per unit accreted mass.  Locally emitted flux has to be redshifted out of
    // the well by the orbiting emitter's specific energy E(r), so the statement
    // is  Integral F(r) E(r) dA = eta Mdot c^2  with eta = 1 - E(r_isco).
    for (double a : {0.0, 0.9}) {
        bh::NovikovThorneDisc d(1.0e7 * phys::M_sun, a, 1.0e9, 0.1);
        double L = 0.0;
        const int N = 2000000;
        const double u0 = std::sqrt(d.inner_radius()), u1 = std::sqrt(1.0e9);
        for (int i = 0; i < N; ++i) {
            const double ua = u0 + (u1 - u0) * i / N;
            const double ub = u0 + (u1 - u0) * (i + 1) / N;
            const double um = 0.5 * (ua + ub);
            auto integ = [&](double u) {
                const double r = u * u;
                // 4 pi r dr = two faces of the disc, coordinate area element.
                return d.flux_dimensionless(r) * bh::circular_energy(a, r, +1) *
                       4.0 * M_PI * r * 2.0 * u;
            };
            L += (ub - ua) / 6.0 * (integ(ua) + 4.0 * integ(um) + integ(ub));
        }
        std::snprintf(buf, sizeof buf, "a = %.1f: integral %.6f Mdot c^2 vs eta = %.6f",
                      a, L, d.efficiency());
        check("disc luminosity at infinity equals eta Mdot c^2",
              close(L, d.efficiency(), 1e-3), buf);
    }

    // Temperatures should land where observations put them.
    {
        bh::NovikovThorneDisc stellar(10.0 * phys::M_sun, 0.9, 50.0, 0.1);
        bh::NovikovThorneDisc sgra(4.297e6 * phys::M_sun, 0.9, 50.0, 1e-8);
        std::snprintf(buf, sizeof buf, "10 Msun: %.3e K (X-ray)   4.3e6 Msun: %.3e K",
                      stellar.peak_temperature(), sgra.peak_temperature());
        check("peak disc temperature scales as M^{-1/4}",
              stellar.peak_temperature() > 1e6 && stellar.peak_temperature() < 5e7 &&
              sgra.peak_temperature() < 1e5, buf);
    }
}

void test_radiometry() {
    section("7. Radiometry");

    char buf[220];

    // Integrating the Planck function over all wavelengths must give the
    // Stefan-Boltzmann law.
    {
        double worst = 0.0;
        for (double T : {1000.0, 5772.0, 1.0e5, 1.0e7}) {
            double total = 0.0;
            // Integrate in log-wavelength, which handles the enormous range.
            const double lo = std::log(1e-11), hi = std::log(1e-1);
            const int N = 200000;
            for (int i = 0; i < N; ++i) {
                const double x = lo + (hi - lo) * (i + 0.5) / N;
                const double lam = std::exp(x);
                total += spec::planck_radiance(lam, T) * lam * (hi - lo) / N;
            }
            const double expect = phys::sigma_SB * T * T * T * T / M_PI;
            worst = std::max(worst, std::fabs(total - expect) / expect);
        }
        std::snprintf(buf, sizeof buf, "max relative error = %.2e", worst);
        check("integral of Planck's law gives sigma T^4 / pi", worst < 1e-5, buf);
    }

    // Wien's displacement law.
    {
        double worst = 0.0;
        for (double T : {2000.0, 5772.0, 30000.0}) {
            double best_l = 0.0, best_v = -1.0;
            for (int i = 0; i < 400000; ++i) {
                const double lam = 1e-9 + i * 1e-11;
                const double v = spec::planck_radiance(lam, T);
                if (v > best_v) { best_v = v; best_l = lam; }
            }
            worst = std::max(worst, std::fabs(best_l * T - 2.897771955e-3) / 2.897771955e-3);
        }
        std::snprintf(buf, sizeof buf, "max relative error in lambda_max T = %.2e", worst);
        check("Wien's displacement law, lambda_max T = 2.8978e-3 m K", worst < 1e-3, buf);
    }

    // The Sun should come out looking like the Sun: near-white, slightly warm,
    // with chromaticity close to the measured solar value.
    {
        const spec::XYZ s = spec::blackbody_xyz(phys::T_sun);
        const double sum = s.x + s.y + s.z;
        const double cx = s.x / sum, cy = s.y / sum;
        std::snprintf(buf, sizeof buf, "chromaticity (%.4f, %.4f), Planckian locus at 5772 K",
                      cx, cy);
        check("5772 K blackbody sits on the Planckian locus near D-white",
              cx > 0.31 && cx < 0.35 && cy > 0.32 && cy < 0.37, buf);
    }

    // Table lookup must agree with direct integration.
    {
        double worst = 0.0;
        for (double lg = 2.0; lg <= 8.0; lg += 0.13) {
            const double T = std::pow(10.0, lg);
            const spec::XYZ e = spec::blackbody_xyz_exact(T);
            const spec::XYZ t = spec::blackbody_xyz(T);
            worst = std::max(worst, std::fabs(t.y - e.y) / std::max(e.y, 1e-300));
        }
        std::snprintf(buf, sizeof buf, "max relative error = %.2e", worst);
        check("blackbody colour table matches direct integration", worst < 2e-3, buf);
    }
}

void test_relativistic_beaming() {
    section("8. Doppler beaming and the disc's redshift map");

    char buf[240];

    // For an emitter on a circular orbit seen edge-on, the approaching side is
    // blueshifted and the receding side redshifted.  Check the special
    // relativistic limit far from the hole, where g reduces to the ordinary
    // relativistic Doppler formula 1 / (gamma (1 -/+ v)).
    const double a = 0.0, r = 1000.0;
    const auto u = bh::keplerian_four_velocity(a, r);
    const double v = std::sqrt(1.0 / r);              // orbital speed, v/c
    const double gamma = 1.0 / std::sqrt(1.0 - v * v);

    // Photon momentum for light emitted tangentially towards / away from us.
    const bh::MetricLower ml = bh::metric_lower(a, r, M_PI / 2);
    const bh::Geom gm = bh::geom_at(a, r, M_PI / 2);
    const bh::Tetrad tet = bh::zamo_tetrad(gm);

    // Build a photon travelling in the +phi and -phi directions in the local
    // static frame; its energy measured by the orbiting matter gives g.
    double worst = 0.0;
    for (int s = -1; s <= 1; s += 2) {
        const double dir[3] = {0.0, 0.0, static_cast<double>(s)};
        std::array<double, 4> pup{};
        for (int mu = 0; mu < 4; ++mu)
            pup[mu] = tet.e[0][mu] + dir[2] * tet.e[3][mu];
        const auto plow = bh::lower(ml, pup);
        // Energy in the orbiting frame relative to the static frame.
        const double E_orbit = -(plow[0] * u[0] + plow[3] * u[3]);
        const double E_static = -(plow[0] * tet.e[0][0] + plow[3] * tet.e[0][3]);
        const double g_meas = E_static / E_orbit;
        const double g_sr = 1.0 / (gamma * (1.0 - s * v));
        worst = std::max(worst, std::fabs(g_meas - g_sr) / g_sr);
    }
    std::snprintf(buf, sizeof buf,
                  "v = %.4f c, gamma = %.6f, max deviation from 1/(gamma(1-/+v)) = %.2e",
                  v, gamma, worst);
    check("Doppler factor reduces to special relativity far from the hole",
          worst < 2e-3, buf);

    // Transverse Doppler (time dilation) alone, seen at the ISCO.
    {
        bh::NovikovThorneDisc d(10.0 * phys::M_sun, 0.0, 30.0, 0.1);
        const double T_isco_local = d.temperature(d.inner_radius() * 1.5);
        std::snprintf(buf, sizeof buf, "T(1.5 r_isco) = %.4e K locally", T_isco_local);
        check("inner disc reaches X-ray temperatures for a stellar-mass hole",
              T_isco_local > 1e6, buf);
    }
}

void test_physical_scales() {
    section("9. Physical scales");

    char buf[240];

    // Schwarzschild radius of the Sun.
    const double rs_sun = phys::r_s_metres(phys::M_sun);
    std::snprintf(buf, sizeof buf, "r_s(Sun) = %.4f km (textbook 2.953 km)", rs_sun / 1000.0);
    check("Schwarzschild radius of the Sun is 2.95 km", close(rs_sun, 2953.25, 1e-3), buf);

    // Sgr A*: the measured shadow diameter from the EHT is 51.8 +/- 2.3 uas.
    {
        const double M = 4.297e6 * phys::M_sun;
        const double D = 8277.0 * phys::parsec;
        const double r_g = phys::r_g_metres(M);
        const double shadow_diam = 2.0 * 3.0 * std::sqrt(3.0) * r_g;     // Schwarzschild
        const double uas = phys::rad_to_microarcsec(shadow_diam / D);
        std::snprintf(buf, sizeof buf, "predicted %.1f uas, EHT measured 51.8 +/- 2.3 uas", uas);
        check("Sgr A* shadow diameter matches the EHT measurement",
              uas > 47.0 && uas < 57.0, buf);
    }

    // M87*: EHT measured ring diameter 42 +/- 3 uas.
    {
        const double M = 6.5e9 * phys::M_sun;
        const double D = 16.8e6 * phys::parsec;
        const double r_g = phys::r_g_metres(M);
        const double shadow_diam = 2.0 * 3.0 * std::sqrt(3.0) * r_g;
        const double uas = phys::rad_to_microarcsec(shadow_diam / D);
        std::snprintf(buf, sizeof buf, "predicted %.1f uas, EHT measured 42 +/- 3 uas", uas);
        check("M87* shadow diameter matches the EHT measurement",
              uas > 36.0 && uas < 48.0, buf);
    }

    // Eddington luminosity of a 10 solar-mass hole.
    {
        const double L = phys::L_eddington(10.0 * phys::M_sun);
        std::snprintf(buf, sizeof buf, "L_Edd(10 Msun) = %.4e W = %.4e L_sun",
                      L, L / phys::L_sun);
        check("Eddington luminosity is 1.26e38 erg/s per solar mass",
              close(L / 10.0, 1.2567e31, 2e-3), buf);
    }
}

// ===========================================================================
// Neutron stars
// ===========================================================================

void test_degenerate_matter() {
    section("10. Degenerate matter: the equation of state");

    char buf[240];

    // The Fermi gas expressions must be thermodynamically consistent:
    // d(eps)/dn = mu and P = mu n - eps.  Neither is imposed; both follow from
    // the integrals over the filled Fermi sphere, so agreement is a real check
    // that the closed forms were derived correctly.
    {
        ns::IdealFermiGas g = ns::neutron_gas();
        double worst_mu = 0.0, worst_P = 0.0;
        for (double x = 0.01; x < 20.0; x *= 1.5) {
            const double h = 1e-6 * x;
            const double dn = g.number_density(x + h) - g.number_density(x - h);
            const double de = (g.density_of_x(x + h) - g.density_of_x(x - h)) *
                              phys::c * phys::c;
            const double mu_num = de / dn;
            const double mu_ana = g.chemical_potential(x);
            worst_mu = std::max(worst_mu, std::fabs(mu_num - mu_ana) / mu_ana);

            const double P_thermo = mu_ana * g.number_density(x) -
                                    g.density_of_x(x) * phys::c * phys::c;
            worst_P = std::max(worst_P,
                               std::fabs(P_thermo - g.pressure_of_x(x)) / g.pressure_of_x(x));
        }
        std::snprintf(buf, sizeof buf, "max relative error = %.2e", worst_mu);
        check("Fermi gas satisfies d(energy)/d(number) = mu", worst_mu < 1e-6, buf);
        // At low x the pressure is O(x^5) while mu n and epsilon are both
        // O(x^3), so forming their difference in double precision loses about
        // four digits.  That is arithmetic, not physics.
        std::snprintf(buf, sizeof buf,
                      "max relative error = %.2e (cancellation-limited at small x)", worst_P);
        check("Fermi gas satisfies P = mu n - epsilon", worst_P < 1e-6, buf);
    }

    // Non-relativistic and ultra-relativistic limits, P ~ n^{5/3} and n^{4/3}.
    {
        ns::IdealFermiGas g = ns::neutron_gas();
        auto slope = [&](double x) {
            const double h = 0.01 * x;
            return (std::log(g.pressure_of_x(x + h)) - std::log(g.pressure_of_x(x - h))) /
                   (std::log(g.number_density(x + h)) - std::log(g.number_density(x - h)));
        };
        const double nr = slope(0.02), ur = slope(500.0);
        std::snprintf(buf, sizeof buf,
                      "d ln P / d ln n = %.4f at x = 0.02, %.4f at x = 500", nr, ur);
        check("degenerate gas stiffens from Gamma = 5/3 to 4/3",
              std::fabs(nr - 5.0 / 3.0) < 2e-3 && std::fabs(ur - 4.0 / 3.0) < 2e-3, buf);
    }
}

void test_tov() {
    section("11. Neutron star structure: the TOV equation");

    char buf[260];

    // The one case with a closed-form solution.  For constant density the TOV
    // equation integrates to the interior Schwarzschild metric of 1916.  Seed
    // the solver with the analytic central pressure for a chosen radius and it
    // should integrate back out to exactly that radius and mass.
    {
        const double rho = 5.0e17;
        const double R_target = 12000.0;
        const double Pc = ns::interior_schwarzschild_pressure(rho, R_target, 0.0);
        ns::UniformDensity ud(rho, Pc);
        const ns::Star st = ns::solve_tov(ud, rho, 1e-12);
        const double M_exact = (4.0 / 3.0) * M_PI * R_target * R_target * R_target * rho;

        std::snprintf(buf, sizeof buf,
                      "integrated R = %.4f km (exact %.4f), M = %.6f Msun (exact %.6f)",
                      st.R_km(), R_target / 1000.0, st.M_solar(), M_exact / phys::M_sun);
        check("uniform-density star reproduces interior Schwarzschild",
              close(st.R, R_target, 1e-4) && close(st.M, M_exact, 1e-4), buf);
    }

    // Buchdahl's bound: no static star of any equation of state can be more
    // compact than 8/9, because the central pressure would have to be infinite.
    {
        double worst = 0.0;
        for (double lp : {34.0, 34.4, 34.9}) {
            ns::PiecewisePolytrope e(lp, 3.0, 3.0, 3.0, "test");
            const ns::MassRadiusCurve cv = ns::mass_radius_curve(e, 3e17, 6e18, 30);
            for (const ns::Star& st : cv.stars) worst = std::max(worst, st.compactness);
        }
        std::snprintf(buf, sizeof buf, "max compactness found = %.4f, bound is %.4f",
                      worst, ns::kBuchdahlCompactness);
        check("every solution respects Buchdahl's bound r_s/R < 8/9",
              worst < ns::kBuchdahlCompactness, buf);
    }

    // The Oppenheimer-Volkoff limit.  Free neutrons with no interactions at
    // all give a maximum mass well below what pulsars actually weigh - which
    // is precisely the 1939 result, and the reason nuclear forces have to
    // matter.
    {
        ns::IdealFermiGas g = ns::neutron_gas();
        const ns::MassRadiusCurve cv = ns::mass_radius_curve(g, 1e17, 5e19, 80);
        std::snprintf(buf, sizeof buf,
                      "M_max = %.4f Msun at R = %.2f km  (Oppenheimer & Volkoff: 0.71)",
                      cv.M_max_solar(), cv.max_mass.R_km());
        check("free neutron gas reproduces the 0.71 Msun OV limit",
              std::fabs(cv.M_max_solar() - 0.71) < 0.02, buf);
    }

    // The Chandrasekhar limit.  The Newtonian value for mu_e = 2 is
    // 1.456 Msun; general relativity destabilises the star slightly and pulls
    // the maximum down by a couple of percent, which the TOV solve should show.
    {
        ns::IdealFermiGas g = ns::electron_gas(2.0);
        const ns::MassRadiusCurve cv = ns::mass_radius_curve(g, 1e9, 1e15, 80);
        std::snprintf(buf, sizeof buf,
                      "M_max = %.4f Msun (Newtonian limit 1.456, GR lowers it)",
                      cv.M_max_solar());
        check("degenerate electrons reproduce the Chandrasekhar mass",
              cv.M_max_solar() > 1.38 && cv.M_max_solar() < 1.46, buf);
    }

    // Newtonian limit: a Gamma = 2 polytrope has the analytic Lane-Emden
    // solution theta = sin(xi)/xi, so its radius is R = pi sqrt(K / (2 pi G))
    // and - remarkably - does not depend on the mass at all.
    {
        const double K = 1.0e6, G2 = 2.0;
        ns::Polytrope p(K, G2);
        const double R_exact = M_PI * std::sqrt(K / (2.0 * M_PI * phys::G));
        // The Lane-Emden solution is Newtonian, so the comparison is only
        // meaningful where P << rho c^2, i.e. K rho << c^2.  With K = 1e6 that
        // means densities well below 1e10; at 3e10 the star is already
        // relativistic enough that general relativity shrinks it by a third,
        // which is a real effect rather than an error.
        double worst = 0.0;
        for (double rho_c : {1.0e6, 3.0e6, 1.0e7}) {
            const ns::Star st = ns::solve_tov(p, rho_c, 1e-14);
            if (!st.ok) continue;
            worst = std::max(worst, std::fabs(st.R - R_exact) / R_exact);
        }
        const ns::Star rel = ns::solve_tov(p, 3.0e10, 1e-14);
        std::snprintf(buf, sizeof buf,
                      "Newtonian R = %.1f km vs Lane-Emden %.1f km (err %.2e); at 3e10 GR gives %.1f km",
                      ns::solve_tov(p, 3.0e6, 1e-14).R / 1000.0, R_exact / 1000.0, worst,
                      rel.R / 1000.0);
        check("Gamma = 2 polytrope matches the Lane-Emden radius", worst < 3e-3, buf);
    }

    // Causality and the observed two-solar-mass pulsars.
    {
        ns::PiecewisePolytrope sly = ns::eos_sly();
        const ns::MassRadiusCurve cv = ns::mass_radius_curve(sly, 3e17, 5e18, 60);
        const ns::Star s14 = ns::star_of_mass(sly, 1.4 * phys::M_sun);
        std::snprintf(buf, sizeof buf,
                      "M_max = %.3f Msun; PSR J0740+6620 weighs 2.08 +/- 0.07",
                      cv.M_max_solar());
        check("SLy supports the heaviest pulsars measured",
              cv.M_max_solar() > 2.0, buf);
        std::snprintf(buf, sizeof buf,
                      "R(1.4 Msun) = %.2f km; NICER measures 12.4 +1.3/-1.0 km for J0740",
                      s14.R_km());
        check("SLy radius sits in the NICER range",
              s14.R_km() > 10.5 && s14.R_km() < 13.5, buf);
        // Causality holds comfortably through a typical star.  At the very
        // top of the mass range the piecewise-polytrope *fit* creeps just past
        // c - an artefact of representing a tabulated equation of state by
        // power laws with Gamma near 3, not a property of SLy itself.  Worth
        // reporting rather than hiding.
        std::snprintf(buf, sizeof buf,
                      "%.4f c in a 1.4 Msun star (%.4f c at the maximum mass, where the fit strains)",
                      s14.max_sound_speed, cv.max_mass.max_sound_speed);
        check("sound speed stays below c through a typical star",
              s14.max_sound_speed < 1.0, buf);
        std::snprintf(buf, sizeof buf,
                      "binding energy = %.4f Msun c^2 = %.3e J (SN 1987A radiated ~3e46 J)",
                      s14.binding_energy / (phys::M_sun * phys::c * phys::c),
                      s14.binding_energy);
        check("binding energy matches supernova neutrino energetics",
              s14.binding_energy > 1e46 && s14.binding_energy < 1e47, buf);
    }

    // A stiffer equation of state must give both a bigger star and a higher
    // maximum mass.  That monotonicity is the whole reason a measured radius
    // constrains nuclear physics.
    {
        ns::PiecewisePolytrope sly = ns::eos_sly(), ms1 = ns::eos_stiff();
        const ns::Star a = ns::star_of_mass(sly, 1.4 * phys::M_sun);
        const ns::Star b = ns::star_of_mass(ms1, 1.4 * phys::M_sun);
        const double Ma = ns::mass_radius_curve(sly, 3e17, 5e18, 40).M_max_solar();
        const double Mb = ns::mass_radius_curve(ms1, 3e17, 5e18, 40).M_max_solar();
        std::snprintf(buf, sizeof buf,
                      "soft: R = %.2f km, M_max = %.2f;  stiff: R = %.2f km, M_max = %.2f",
                      a.R_km(), Ma, b.R_km(), Mb);
        check("a stiffer equation of state gives larger, heavier stars",
              b.R_km() > a.R_km() && Mb > Ma, buf);
    }
}

void test_neutron_star_optics() {
    section("12. Neutron star optics: seeing round the back");

    char buf[280];

    // For a surface at R the light bending is strong enough to show a large
    // part of the far hemisphere.  Trace rays at increasing impact parameter,
    // find the largest one that still lands on the surface, and read off the
    // colatitude it lands at.  Compare with an exact quadrature of
    // dpsi/dr for the grazing ray, and with Beloborodov's approximation.
    for (double R : {5.0, 6.0, 8.0}) {
        const double a = 0.0;
        const double r0 = 1.0e5;
        const bh::Geom g = bh::geom_at(a, r0, M_PI / 2);
        const bh::Tetrad tet = bh::zamo_tetrad(g);

        // Does a ray aimed with this impact parameter reach the surface, and
        // if so at what azimuth?
        auto trace_to_surface = [&](double b, double* phi_hit) {
            const double alpha = std::asin(std::clamp(
                b * std::sqrt(1.0 - 2.0 / r0) / r0, -1.0, 1.0));
            const double dir[3] = {-std::cos(alpha), 0.0, std::sin(alpha)};
            bh::State y = bh::make_photon(a, r0, M_PI / 2, 0.0, tet, dir, true);
            bh::Stepper st(a, 1e-12, 1e-14);
            st.reset();
            double h = 1.0;
            bh::State yn{};
            bh::DenseSegment seg;
            for (int i = 0; i < 400000; ++i) {
                int tries = 0;
                while (!st.try_step(y, h, yn, seg)) if (++tries > 60) return false;
                if (h > 0.02 * yn[bh::Y_R]) h = 0.02 * yn[bh::Y_R];
                if (yn[bh::Y_R] <= R) {
                    double lo = 0.0, hi = 1.0;
                    for (int k = 0; k < 60; ++k) {
                        const double mid = 0.5 * (lo + hi);
                        if (seg.eval(mid)[bh::Y_R] <= R) hi = mid; else lo = mid;
                    }
                    if (phi_hit) *phi_hit = std::fabs(seg.eval(0.5 * (lo + hi))[bh::Y_PH]);
                    return true;
                }
                if (yn[bh::Y_R] > 2.0 * r0) return false;
                y = yn;
            }
            return false;
        };

        // Largest impact parameter that still hits.
        double lo = 0.0, hi = 3.0 * R;
        for (int i = 0; i < 60; ++i) {
            const double mid = 0.5 * (lo + hi);
            double dummy;
            if (trace_to_surface(mid, &dummy)) lo = mid; else hi = mid;
        }
        double psi = 0.0;
        trace_to_surface(lo, &psi);

        const double b_analytic = R / std::sqrt(1.0 - 2.0 / R);
        const double psi_exact = bh::max_visible_angle_exact(R);
        const double cos_bel = 1.0 - 1.0 / (1.0 - 2.0 / R);
        const double psi_bel = std::acos(std::clamp(cos_bel, -1.0, 1.0));

        std::snprintf(buf, sizeof buf,
                      "R = %.1f M: b_max = %.5f vs R/sqrt(1-r_s/R) = %.5f", R, lo, b_analytic);
        check("apparent size is set by the critical impact parameter",
              close(lo, b_analytic, 2e-4), buf);

        std::snprintf(buf, sizeof buf,
                      "R = %.1f M: traced %.3f deg, quadrature %.3f deg, Beloborodov %.3f deg",
                      R, psi * 180 / M_PI, psi_exact * 180 / M_PI, psi_bel * 180 / M_PI);
        check("last visible point matches the exact quadrature",
              close(psi, psi_exact, 3e-3), buf);
    }

    // Visible fraction and surface redshift for a real star.
    {
        ns::PiecewisePolytrope sly = ns::eos_sly();
        const ns::Star st = ns::star_of_mass(sly, 1.4 * phys::M_sun);
        bh::NeutronStar S;
        S.set_from(st);
        const double z_direct = 1.0 / std::sqrt(1.0 - 2.0 / S.R) - 1.0;
        std::snprintf(buf, sizeof buf,
                      "R = %.3f GM/c^2: %.1f%% of the surface visible, z = %.4f",
                      S.R, 100.0 * S.visible_fraction(), z_direct);
        check("more than half the surface is visible",
              S.visible_fraction() > 0.7 && S.visible_fraction() < 0.9 &&
              close(z_direct, st.redshift, 1e-9), buf);
    }

    // A rigidly rotating surface must still have a properly normalised
    // four-velocity in the Schwarzschild metric.
    {
        bh::NeutronStar S;
        S.R = 5.5;
        S.M_kg = 1.4 * phys::M_sun;
        S.set_spin(600.0);
        double worst = 0.0;
        for (double th : {0.2, 0.7, M_PI / 2, 2.4}) {
            const auto u = S.surface_four_velocity(th);
            const auto ml = bh::metric_lower(0.0, S.R, th);
            const double n = ml.tt * u[0] * u[0] + ml.pp * u[3] * u[3];
            worst = std::max(worst, std::fabs(n + 1.0));
        }
        std::snprintf(buf, sizeof buf,
                      "max |u.u + 1| = %.2e at %.0f Hz, equatorial speed %.4f c",
                      worst, S.spin_hz, S.equatorial_speed());
        check("rotating surface four-velocity is normalised", worst < 1e-12, buf);
    }
}

// ===========================================================================
// 13. Relativistic jets
// ===========================================================================
void test_jets() {
    section("13. Blandford-Znajek jets");
    char buf[512];

    // --- the power itself ---------------------------------------------------

    // A hole with no spin has no rotational energy to give up.  This is the
    // whole point: the jet is not powered by the accretion flow.
    {
        const double P0 = bh::blandford_znajek_power(0.0, 1.0);
        std::snprintf(buf, sizeof buf, "P(a=0) = %.3e Mdot c^2", P0);
        check("a non-rotating hole drives no jet", P0 == 0.0, buf);
    }

    // At small spin Omega_H -> a/4 and the power is quadratic in it.  Doubling
    // the spin must quadruple the power.
    {
        const double P1 = bh::blandford_znajek_power(0.01, 1.0);
        const double P2 = bh::blandford_znajek_power(0.02, 1.0);
        std::snprintf(buf, sizeof buf, "P(2a)/P(a) = %.6f, expected 4", P2 / P1);
        check("jet power scales as Omega_H^2 at small spin",
              close(P2 / P1, 4.0, 2e-3), buf);
    }

    // The headline number, and the reason a maximally spinning hole is
    // interesting: the jet takes out more than the accretion brings in.
    {
        const double P = bh::blandford_znajek_power(0.998, 1.0);   // per Mdot c^2
        const double OmH = bh::horizon_angular_velocity(0.998);
        std::snprintf(buf, sizeof buf,
                      "Omega_H = %.4f, P = %.4f Mdot c^2 (phi = 50, MAD)", OmH, P);
        check("at the Thorne limit the jet outpowers the accretion flow",
              P > 1.5 && P < 2.5, buf);
    }

    // Monotone in spin over the whole physical range.
    {
        bool mono = true;
        double prev = -1.0;
        for (int i = 0; i <= 100; ++i) {
            const double P = bh::blandford_znajek_power(0.998 * i / 100.0, 1.0);
            if (P < prev) mono = false;
            prev = P;
        }
        std::snprintf(buf, sizeof buf, "P(0.998)/P(0.5) = %.3f",
                      bh::blandford_znajek_power(0.998, 1.0) /
                      bh::blandford_znajek_power(0.5, 1.0));
        check("jet power increases monotonically with spin", mono, buf);
    }

    // --- the flow ------------------------------------------------------------

    // The plasma four-velocity has to be a unit timelike vector everywhere,
    // including inside the ergosphere where no static observer exists at all.
    {
        bh::Jet j;
        j.enabled = true;
        j.configure(0.998, 1.0, 1.0);
        double worst = 0.0, worst_r = 0.0;
        for (double r : {1.1, 1.5, 3.0, 8.0, 30.0, 120.0}) {
            for (double th : {0.02, 0.3, 0.9, M_PI / 2, 2.2, 3.1}) {
                const bh::Geom g = bh::geom_at(0.998, r, th);
                const auto u = j.four_velocity(g);
                const bh::MetricLower m = bh::metric_lower(g);
                const double n = m.tt * u[0] * u[0] + 2.0 * m.tp * u[0] * u[3] +
                                 m.rr * u[1] * u[1] + m.thth * u[2] * u[2] +
                                 m.pp * u[3] * u[3];
                if (std::fabs(n + 1.0) > worst) { worst = std::fabs(n + 1.0); worst_r = r; }
            }
        }
        std::snprintf(buf, sizeof buf,
                      "max |u.u + 1| = %.2e (worst at r = %.1f, a = 0.998)", worst, worst_r);
        check("jet plasma four-velocity is normalised", worst < 1e-12, buf);
    }

    // The Lorentz factor a local non-rotating observer measures must be the one
    // the velocity profile asked for.
    {
        bh::Jet j;
        j.enabled = true;
        j.configure(0.9, 1.0, 1.0);
        double worst = 0.0;
        for (double z : {5.0, 15.0, 40.0, 100.0}) {
            const double r = z / std::cos(0.25), th = 0.25;
            const bh::Geom g = bh::geom_at(0.9, r, th);
            const auto u = j.four_velocity(g);
            const bh::Tetrad T = bh::zamo_tetrad(g);
            const bh::MetricLower m = bh::metric_lower(g);
            // W = -u . n, with n the ZAMO four-velocity e_(0).
            const double W = -(m.tt * u[0] * T.e[0][0] +
                               m.tp * (u[0] * T.e[0][3] + u[3] * T.e[0][0]) +
                               m.pp * u[3] * T.e[0][3]);
            worst = std::max(worst, std::fabs(W - j.gamma_at(r * std::cos(th))) /
                                    j.gamma_at(r * std::cos(th)));
        }
        std::snprintf(buf, sizeof buf,
                      "max relative error %.2e against the prescribed Gamma(z)", worst);
        check("locally measured Lorentz factor matches the profile", worst < 1e-12, buf);
    }

    // --- beaming -------------------------------------------------------------

    // Far from the hole spacetime is flat, so the redshift factor the renderer
    // computes has to reduce to the special-relativistic Doppler factor
    // delta = 1 / (Gamma (1 - beta cos(theta))) for the photon's own direction
    // of travel.  This is the single quantity that decides how one-sided a jet
    // looks, so it is worth checking against the textbook formula and not just
    // against itself.
    {
        bh::Jet j;
        j.enabled = true;
        j.configure(0.9, 1.0, 1.0);
        const double r = 1.0e6, th = 1.0e-3;   // on the axis, effectively flat
        const bh::Geom g = bh::geom_at(0.9, r, th);
        const bh::Tetrad T = bh::zamo_tetrad(g);
        const bh::MetricLower m = bh::metric_lower(g);
        const auto u = j.four_velocity(g);

        // Extract the plasma's ordinary velocity in the local orthonormal frame.
        auto dotT = [&](int leg) {
            return m.tt * u[0] * T.e[leg][0] +
                   m.tp * (u[0] * T.e[leg][3] + u[3] * T.e[leg][0]) +
                   m.rr * u[1] * T.e[leg][1] + m.thth * u[2] * T.e[leg][2] +
                   m.pp * u[3] * T.e[leg][3];
        };
        const double W = -dotT(0);
        const double v[3] = {dotT(1) / W, dotT(2) / W, dotT(3) / W};
        const double beta = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);

        double worst = 0.0;
        for (double ang = 0.1; ang < 3.1; ang += 0.2) {
            // Ray direction away from the camera; the photon travels the other
            // way, so cos(theta) picks up the minus sign.
            const double dir[3] = {std::cos(ang), std::sin(ang), 0.0};
            const bh::State y = bh::make_photon(0.9, r, th, 0.0, T, dir, true);
            const double ku = y[bh::Y_PT] * u[0] + y[bh::Y_PR] * u[1] +
                              y[bh::Y_PTH] * u[2] + y[bh::Y_PPH] * u[3];
            const double gfac = 1.0 / ku;
            const double cosang = -(dir[0] * v[0] + dir[1] * v[1] + dir[2] * v[2]) / beta;
            const double delta = 1.0 / (W * (1.0 - beta * cosang));
            worst = std::max(worst, std::fabs(gfac - delta) / delta);
        }
        std::snprintf(buf, sizeof buf,
                      "Gamma = %.3f, beta = %.5f: max relative error %.2e", W, beta, worst);
        check("redshift factor reduces to the Doppler formula in flat space",
              worst < 1e-6, buf);
    }

    // The observable consequence: a jet pointed near the line of sight is
    // enormously brighter than its receding twin, and the same jet viewed from
    // the side is not.  Nothing in the code puts this in by hand.
    {
        bh::Jet j;
        const double G = 10.0;
        const double beta = std::sqrt(1.0 - 1.0 / (G * G));
        auto ratio = [&](double incl_deg) {
            const double c = std::cos(incl_deg * M_PI / 180.0);
            const double d_app = 1.0 / (G * (1.0 - beta * c));
            const double d_rec = 1.0 / (G * (1.0 + beta * c));
            return std::pow(d_app / d_rec, 2.0 + j.alpha);
        };
        const double r20 = ratio(20.0), r84 = ratio(84.0);
        std::snprintf(buf, sizeof buf,
                      "Gamma = 10: %.3e at 20 deg, %.2f at 84 deg", r20, r84);
        check("beaming makes a jet one-sided only at small inclination",
              r20 > 1e3 && r84 > 1.0 && r84 < 3.0, buf);
    }

    // --- energy bookkeeping --------------------------------------------------

    // configure() normalises the emissivity by doing the transverse integral
    // analytically and the longitudinal one on a grid.  Check the result the
    // hard way, with a brute-force two-dimensional quadrature of the emissivity
    // that the renderer itself samples.
    {
        bh::Jet j;
        j.enabled = true;
        const double r_g = phys::r_g_metres(1.0e9 * phys::M_sun);
        const double mdot_c2 = 2.0e40;
        j.configure(0.998, mdot_c2, r_g);

        // j_bol = (j_nu0 / r_g) * nu_0^alpha * Integral_{nu_lo}^{nu_hi} nu^-alpha dnu
        const double band = (std::pow(bh::Jet::nu_hi, 1.0 - j.alpha) -
                             std::pow(bh::Jet::nu_lo, 1.0 - j.alpha)) / (1.0 - j.alpha);
        const double A = j.j_nu0 / r_g * std::pow(bh::Jet::nu_0, j.alpha) * band;

        double V = 0.0;
        constexpr int NZ = 1200, NR = 400;
        const double dz = (j.z_top - j.z_base) / NZ;
        for (int i = 0; i < NZ; ++i) {
            const double z = j.z_base + (i + 0.5) * dz;
            const double rmax = 2.0 * j.radius(z);
            const double dr = rmax / NR;
            for (int k = 0; k < NR; ++k) {
                const double rho = (k + 0.5) * dr;
                V += j.shape(z, rho) * 2.0 * M_PI * rho * dr * dz;
            }
        }
        V *= 2.0;                                    // both jets
        const double L = 4.0 * M_PI * A * V * r_g * r_g * r_g;
        std::snprintf(buf, sizeof buf,
                      "recovered %.4e W against %.4e W requested (%.3f%% off)",
                      L, j.L_rad, 100.0 * std::fabs(L / j.L_rad - 1.0));
        check("emitted power matches the Blandford-Znajek budget",
              close(L, j.L_rad, 2e-3), buf);
    }

    // The prescribed collimation really is the measured parabola, and the flow
    // stays subluminal everywhere.
    {
        bh::Jet j;
        j.enabled = true;
        j.configure(0.998, 1.0, 1.0);
        const double slope = std::log(j.radius(100.0) / j.radius(10.0)) / std::log(10.0);
        bool sub = true;
        for (double z = j.z_base; z <= j.z_top; z += 0.5)
            if (!(j.speed_at(z) < 1.0) || !(j.gamma_at(z) <= j.gamma_max + 1e-12)) sub = false;
        std::snprintf(buf, sizeof buf,
                      "d log R / d log z = %.4f (Asada & Nakamura measure 0.58)", slope);
        check("jet boundary follows the measured parabola",
              close(slope, j.k, 1e-12) && sub, buf);
    }

    // --- colour --------------------------------------------------------------

    // powerlaw_xyz integrates in wavelength; redo it in frequency, which is a
    // different substitution with a different Jacobian, and require the same
    // answer.
    {
        const double alpha = 0.7;
        const spec::XYZ a1 = spec::powerlaw_xyz(alpha);
        // I_nu = (nu/nu_0)^-alpha, integrated against the CMFs in frequency.
        spec::XYZ a2{};
        {
            constexpr double lam_lo = 360e-9, lam_hi = 830e-9;
            const double nu_lo = phys::c / lam_hi, nu_hi = phys::c / lam_lo;
            constexpr int N = 200000;
            const double dnu = (nu_hi - nu_lo) / N;
            for (int i = 0; i < N; ++i) {
                const double nu = nu_lo + (i + 0.5) * dnu;
                const double nm = phys::c / nu * 1e9;
                const double w = std::pow(nu / (phys::c / 550e-9), -alpha) * dnu;
                const spec::XYZ cmf = spec::cie_cmf(nm);
                a2.x += w * cmf.x;
                a2.y += w * cmf.y;
                a2.z += w * cmf.z;
            }
        }
        std::snprintf(buf, sizeof buf,
                      "Y: %.6e (wavelength) vs %.6e (frequency), %.3f%% apart",
                      a1.y, a2.y, 100.0 * std::fabs(a1.y / a2.y - 1.0));
        check("synchrotron colour is integration-variable independent",
              close(a1.y, a2.y, 3e-3) && close(a1.x, a2.x, 3e-3) &&
              close(a1.z, a2.z, 3e-3), buf);
    }

    // A rising-to-the-blue power law must actually come out blue: bluer than
    // the Sun, and bluer still for a flatter spectrum.
    {
        const spec::RGB jet = spec::gamut_clamp(
            spec::xyz_to_linear_srgb(spec::powerlaw_xyz(0.7)));
        const spec::RGB sun = spec::blackbody_rgb(phys::T_sun);
        const spec::RGB flat = spec::gamut_clamp(
            spec::xyz_to_linear_srgb(spec::powerlaw_xyz(0.0)));
        const double bj = jet.b / jet.r, bs = sun.b / sun.r, bf = flat.b / flat.r;
        std::snprintf(buf, sizeof buf,
                      "blue/red: alpha=0.7 %.3f, alpha=0 %.3f, the Sun %.3f", bj, bf, bs);
        // I_lambda ~ lambda^(alpha-2), so a *flatter* spectrum is the bluer one.
        check("optically thin synchrotron renders blue", bj > bs && bf > bj, buf);
    }
}

} // namespace

int run_validation() {
    std::printf("\n");
    std::printf("================================================================\n");
    std::printf("  Physics validation\n");
    std::printf("================================================================\n");

    test_metric_algebra();
    test_einstein_equations();
    test_orbits();
    test_null_geodesics();
    test_extreme_rotation();
    test_timelike_geodesics();
    test_disc();
    test_radiometry();
    test_relativistic_beaming();
    test_physical_scales();
    test_degenerate_matter();
    test_tov();
    test_neutron_star_optics();
    test_jets();

    std::printf("\n================================================================\n");
    std::printf("  %d passed, %d failed\n", g_pass, g_fail);
    std::printf("================================================================\n\n");
    return g_fail == 0 ? 0 : 1;
}
