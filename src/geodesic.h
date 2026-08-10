// geodesic.h - Null and timelike geodesics of the Kerr spacetime.
//
// Rather than integrating the second-order geodesic equation
//
//     d2x^mu/dl^2 + Gamma^mu_{ab} dx^a/dl dx^b/dl = 0
//
// (which would need all 40 Christoffel symbols at every step) we integrate the
// equivalent first-order Hamiltonian system.  With the super-Hamiltonian
//
//     H(x, p) = 1/2 g^{mu nu}(x) p_mu p_nu
//
// Hamilton's equations are
//
//     dx^mu / dl = + dH/dp_mu = g^{mu nu} p_nu
//     dp_mu / dl = - dH/dx^mu = -1/2 (d_mu g^{ab}) p_a p_b
//
// which is exactly the geodesic equation written with the covariant momentum.
// Two advantages matter here:
//
//   1. Kerr is stationary and axisymmetric, so d_t g^{ab} = d_phi g^{ab} = 0.
//      The equations for p_t and p_phi therefore read dp_t/dl = dp_phi/dl = 0
//      and the energy E = -p_t and axial angular momentum L_z = p_phi are
//      conserved to machine precision, for free, with no drift at all.
//   2. Only d_r g^{ab} and d_th g^{ab} are ever needed, and kerr.h supplies
//      both analytically.
//
// The remaining two integrals - the norm g^{ab} p_a p_b (= 0 for light) and
// Carter's constant Q - are *not* imposed; they are left free and are used as
// independent accuracy diagnostics.
//
// The integrator is the Dormand-Prince 5(4) embedded Runge-Kutta pair with
// PI step-size control and Hairer's 5th-order dense output, which lets us
// locate events (disc crossings, surface hits) to sub-step accuracy.
#pragma once

#include "kerr.h"

#include <array>
#include <cmath>
#include <algorithm>

namespace bh {

// State vector y = (t, r, theta, phi, p_t, p_r, p_theta, p_phi).
using State = std::array<double, 8>;

enum : int { Y_T = 0, Y_R = 1, Y_TH = 2, Y_PH = 3, Y_PT = 4, Y_PR = 5, Y_PTH = 6, Y_PPH = 7 };

// Right-hand side of Hamilton's equations.
inline State geodesic_rhs(double a, const State& y) {
    const Geom g = geom_at(a, y[Y_R], y[Y_TH]);
    const MetricUpperDeriv m = metric_upper_deriv(g);

    const double pt = y[Y_PT], pr = y[Y_PR], pth = y[Y_PTH], pph = y[Y_PPH];

    State d{};
    // dx^mu/dl = g^{mu nu} p_nu
    d[Y_T]  = m.g.tt * pt + m.g.tp * pph;
    d[Y_R]  = m.g.rr * pr;
    d[Y_TH] = m.g.thth * pth;
    d[Y_PH] = m.g.tp * pt + m.g.pp * pph;

    // dp_mu/dl = -1/2 (d_mu g^{ab}) p_a p_b.  Zero for mu = t, phi.
    d[Y_PT]  = 0.0;
    d[Y_PPH] = 0.0;
    d[Y_PR]  = -0.5 * (m.dr.tt * pt * pt + 2.0 * m.dr.tp * pt * pph +
                       m.dr.rr * pr * pr + m.dr.thth * pth * pth + m.dr.pp * pph * pph);
    d[Y_PTH] = -0.5 * (m.dth.tt * pt * pt + 2.0 * m.dth.tp * pt * pph +
                       m.dth.rr * pr * pr + m.dth.thth * pth * pth + m.dth.pp * pph * pph);
    return d;
}

// ---------------------------------------------------------------------------
// Dormand-Prince 5(4), FSAL, with dense output (Hairer, Norsett & Wanner,
// "Solving Ordinary Differential Equations I", section II.4).
// ---------------------------------------------------------------------------
namespace dp {
inline constexpr double c2 = 1.0 / 5.0, c3 = 3.0 / 10.0, c4 = 4.0 / 5.0, c5 = 8.0 / 9.0;

inline constexpr double a21 = 1.0 / 5.0;
inline constexpr double a31 = 3.0 / 40.0,        a32 = 9.0 / 40.0;
inline constexpr double a41 = 44.0 / 45.0,       a42 = -56.0 / 15.0,      a43 = 32.0 / 9.0;
inline constexpr double a51 = 19372.0 / 6561.0,  a52 = -25360.0 / 2187.0,
                        a53 = 64448.0 / 6561.0,  a54 = -212.0 / 729.0;
inline constexpr double a61 = 9017.0 / 3168.0,   a62 = -355.0 / 33.0,
                        a63 = 46732.0 / 5247.0,  a64 = 49.0 / 176.0,
                        a65 = -5103.0 / 18656.0;
// 5th-order weights (also the 7th stage row: FSAL).
inline constexpr double b1 = 35.0 / 384.0,       b3 = 500.0 / 1113.0,
                        b4 = 125.0 / 192.0,      b5 = -2187.0 / 6784.0,
                        b6 = 11.0 / 84.0;
// Difference between the 5th- and 4th-order weights, for the error estimate.
inline constexpr double e1 = b1 - 5179.0 / 57600.0;
inline constexpr double e3 = b3 - 7571.0 / 16695.0;
inline constexpr double e4 = b4 - 393.0 / 640.0;
inline constexpr double e5 = b5 - (-92097.0 / 339200.0);
inline constexpr double e6 = b6 - 187.0 / 2100.0;
inline constexpr double e7 = 0.0 - 1.0 / 40.0;
// Dense-output coefficients.
inline constexpr double d1 = -12715105075.0 / 11282082432.0;
inline constexpr double d3 = 87487479700.0 / 32700410799.0;
inline constexpr double d4 = -10690763975.0 / 1880347072.0;
inline constexpr double d5 = 701980252875.0 / 199316789632.0;
inline constexpr double d6 = -1453857185.0 / 822651844.0;
inline constexpr double d7 = 69997945.0 / 29380423.0;
} // namespace dp

// Everything needed to evaluate the solution anywhere inside the accepted step.
struct DenseSegment {
    State r1, r2, r3, r4, r5;
    double h = 0.0;

    // Interpolate at fractional position s in [0, 1] across the step.
    State eval(double s) const {
        const double s1 = 1.0 - s;
        State y{};
        for (int i = 0; i < 8; ++i)
            y[i] = r1[i] + s * (r2[i] + s1 * (r3[i] + s * (r4[i] + s1 * r5[i])));
        return y;
    }
};

class Stepper {
public:
    Stepper(double spin, double rtol, double atol)
        : a_(spin), rtol_(rtol), atol_(atol) {}

    // Attempt one step of size h from y. Returns true if accepted; h is updated
    // to the suggested next step size in either case.
    bool try_step(const State& y, double& h, State& y_out, DenseSegment& dense) {
        using namespace dp;
        if (!have_k1_) { k1_ = geodesic_rhs(a_, y); have_k1_ = true; }

        State tmp{};
        for (int i = 0; i < 8; ++i) tmp[i] = y[i] + h * a21 * k1_[i];
        State k2 = geodesic_rhs(a_, tmp);
        for (int i = 0; i < 8; ++i) tmp[i] = y[i] + h * (a31 * k1_[i] + a32 * k2[i]);
        State k3 = geodesic_rhs(a_, tmp);
        for (int i = 0; i < 8; ++i)
            tmp[i] = y[i] + h * (a41 * k1_[i] + a42 * k2[i] + a43 * k3[i]);
        State k4 = geodesic_rhs(a_, tmp);
        for (int i = 0; i < 8; ++i)
            tmp[i] = y[i] + h * (a51 * k1_[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
        State k5 = geodesic_rhs(a_, tmp);
        for (int i = 0; i < 8; ++i)
            tmp[i] = y[i] + h * (a61 * k1_[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] +
                                 a65 * k5[i]);
        State k6 = geodesic_rhs(a_, tmp);

        State ynew{};
        for (int i = 0; i < 8; ++i)
            ynew[i] = y[i] + h * (b1 * k1_[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] +
                                  b6 * k6[i]);
        State k7 = geodesic_rhs(a_, ynew);

        // Scaled error norm (RMS), Hairer's criterion: accept if err <= 1.
        double err2 = 0.0;
        for (int i = 0; i < 8; ++i) {
            const double e = h * (e1 * k1_[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] +
                                  e6 * k6[i] + e7 * k7[i]);
            const double sc = atol_ + rtol_ * std::max(std::fabs(y[i]), std::fabs(ynew[i]));
            const double q = e / sc;
            err2 += q * q;
        }
        double err = std::sqrt(err2 / 8.0);
        if (!std::isfinite(err)) err = 1e10;

        // PI controller (Gustafsson): smoother step sequences than the plain
        // elementary controller, fewer rejections.
        const double alpha = 0.7 / 5.0, beta = 0.4 / 5.0, safety = 0.9;
        double fac = safety * std::pow(std::max(err, 1e-16), -alpha) *
                     std::pow(std::max(err_prev_, 1e-4), beta);
        fac = std::clamp(fac, 0.2, 5.0);

        if (err <= 1.0) {
            // Build the dense-output polynomial before k1 is overwritten.
            dense.h = h;
            for (int i = 0; i < 8; ++i) {
                dense.r1[i] = y[i];
                dense.r2[i] = ynew[i] - y[i];
                dense.r3[i] = h * k1_[i] - dense.r2[i];
                dense.r4[i] = dense.r2[i] - h * k7[i] - dense.r3[i];
                dense.r5[i] = h * (dp::d1 * k1_[i] + dp::d3 * k3[i] + dp::d4 * k4[i] +
                                   dp::d5 * k5[i] + dp::d6 * k6[i] + dp::d7 * k7[i]);
            }
            y_out = ynew;
            k1_ = k7;              // FSAL: last stage becomes the next first stage
            err_prev_ = std::max(err, 1e-4);
            h *= fac;
            return true;
        }
        h *= std::min(fac, 1.0);
        return false;
    }

    // Call after any discontinuous change of state (or a fresh ray).
    void reset() { have_k1_ = false; err_prev_ = 1e-4; }

private:
    double a_, rtol_, atol_;
    State k1_{};
    bool have_k1_ = false;
    double err_prev_ = 1e-4;
};

// ---------------------------------------------------------------------------
// Helpers for setting up and checking geodesics.
// ---------------------------------------------------------------------------

// Given the spatial part of a null momentum in a local orthonormal frame,
// solve g^{ab} p_a p_b = 0 exactly by construction: build p^mu from the tetrad
// and lower it.  `dir` is a unit 3-vector in the frame's (r, th, phi) basis;
// `past_directed` flips the ray so that increasing affine parameter runs
// backwards in coordinate time, which is what backwards ray tracing needs.
inline State make_photon(double a, double r, double th, double phi,
                         const Tetrad& tet, const double dir[3],
                         bool past_directed = true) {
    const Geom g = geom_at(a, r, th);
    const double s = past_directed ? -1.0 : +1.0;

    // p^mu = E (e_(0)^mu + s * dir^i e_(i)^mu), with E = 1 in the local frame.
    std::array<double, 4> pup{};
    for (int mu = 0; mu < 4; ++mu) {
        pup[mu] = tet.e[0][mu] + s * (dir[0] * tet.e[1][mu] + dir[1] * tet.e[2][mu] +
                                      dir[2] * tet.e[3][mu]);
    }
    if (past_directed) for (int mu = 0; mu < 4; ++mu) pup[mu] = -pup[mu];

    const std::array<double, 4> plow = lower(metric_lower(g), pup);
    return {0.0, r, th, phi, plow[0], plow[1], plow[2], plow[3]};
}

// Conserved quantities of a state, for diagnostics.
struct Invariants {
    double E, Lz, Q;
    double norm;      // g^{ab} p_a p_b, exactly 0 for light
    double norm_rel;  // the same, divided by the magnitude of its own terms
};

inline Invariants invariants(double a, const State& y, double m2 = 0.0) {
    Invariants iv;
    iv.E  = -y[Y_PT];
    iv.Lz = y[Y_PPH];
    iv.Q  = carter_constant(a, y[Y_TH], y[Y_PTH], iv.E, iv.Lz, m2);
    const MetricUpper m = metric_upper(a, y[Y_R], y[Y_TH]);
    const std::array<double, 4> p = {y[Y_PT], y[Y_PR], y[Y_PTH], y[Y_PPH]};
    iv.norm = hamiltonian(m, p);
    const double sc = hamiltonian_scale(m, p);
    iv.norm_rel = (sc > 0.0) ? iv.norm / sc : iv.norm;
    return iv;
}

// Pseudo-Cartesian (Kerr-Schild style) embedding of Boyer-Lindquist
// coordinates.  Reduces to ordinary spherical coordinates as a -> 0 and is
// the natural chart for placing distant objects such as the Sun.
inline void bl_to_cartesian(double a, double r, double th, double phi, double out[3]) {
    const double rr = std::sqrt(r * r + a * a);
    const double sn = std::sin(th);
    out[0] = rr * sn * std::cos(phi);
    out[1] = rr * sn * std::sin(phi);
    out[2] = r * std::cos(th);
}

} // namespace bh
