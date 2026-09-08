// kerr.h - The Kerr spacetime in Boyer-Lindquist coordinates.
//
// The Kerr metric is the unique stationary, axisymmetric, asymptotically flat
// vacuum solution of Einstein's field equations
//
//     R_{mu nu} - 1/2 R g_{mu nu} + Lambda g_{mu nu} = 8 pi T_{mu nu}
//
// with Lambda = 0 and T_{mu nu} = 0, i.e.  R_{mu nu} = 0.  (validate.cpp
// verifies R_{mu nu} = 0 numerically straight from the metric components
// defined here, so nothing about the curvature is taken on faith.)
//
// Units: G = c = M = 1.  Lengths are in gravitational radii r_g = GM/c^2.
// Spin is the dimensionless Kerr parameter a = J c / (G M^2), |a| <= 1.
//
// Line element (Boyer-Lindquist, signature -+++):
//
//   ds^2 = -(1 - 2Mr/S) dt^2 - (4 M a r sin^2(th) / S) dt dphi
//          + (S/D) dr^2 + S dth^2 + (A sin^2(th)/S) dphi^2
//
//   S (Sigma) = r^2 + a^2 cos^2(th)
//   D (Delta) = r^2 - 2Mr + a^2
//   A         = (r^2 + a^2)^2 - a^2 D sin^2(th)
#pragma once

#include <array>
#include <cmath>
#include <algorithm>

namespace bh {

// Boyer-Lindquist coordinates are singular on the polar axis (sin(th) = 0).
// Rays that pass exactly through the axis are measure-zero, but finite
// arithmetic still has to be protected; we clamp sin(th) away from zero.
inline constexpr double kSinThetaMin = 1e-9;

// Repeated sub-expressions of the metric, evaluated once per point.
struct Geom {
    double r, th;
    double sn, cs, sn2, cs2;
    double Sigma, Delta, A;
    double a, a2;
};

inline Geom geom_at(double a, double r, double th) {
    Geom g;
    g.a  = a;
    g.a2 = a * a;
    g.r  = r;
    g.th = th;
    g.sn = std::sin(th);
    if (std::fabs(g.sn) < kSinThetaMin) g.sn = (g.sn < 0.0 ? -kSinThetaMin : kSinThetaMin);
    g.cs  = std::cos(th);
    g.sn2 = g.sn * g.sn;
    g.cs2 = g.cs * g.cs;
    g.Sigma = r * r + g.a2 * g.cs2;
    g.Delta = r * r - 2.0 * r + g.a2;
    g.A     = (r * r + g.a2) * (r * r + g.a2) - g.a2 * g.Delta * g.sn2;
    return g;
}

// ---------------------------------------------------------------------------
// Covariant metric g_{mu nu}.  Index order (t, r, th, phi).
// Only the five independent non-zero components are stored.
// ---------------------------------------------------------------------------
struct MetricLower {
    double tt, tp, rr, thth, pp;
};

inline MetricLower metric_lower(const Geom& g) {
    MetricLower m;
    m.tt   = -(1.0 - 2.0 * g.r / g.Sigma);
    m.tp   = -2.0 * g.a * g.r * g.sn2 / g.Sigma;
    m.rr   = g.Sigma / g.Delta;
    m.thth = g.Sigma;
    m.pp   = g.A * g.sn2 / g.Sigma;
    return m;
}

inline MetricLower metric_lower(double a, double r, double th) {
    return metric_lower(geom_at(a, r, th));
}

// Full 4x4 covariant metric (used by the curvature validator).
inline void metric_lower_full(double a, double r, double th, double G[4][4]) {
    MetricLower m = metric_lower(a, r, th);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) G[i][j] = 0.0;
    G[0][0] = m.tt;
    G[1][1] = m.rr;
    G[2][2] = m.thth;
    G[3][3] = m.pp;
    G[0][3] = G[3][0] = m.tp;
}

// ---------------------------------------------------------------------------
// Contravariant metric g^{mu nu} and its analytic partial derivatives.
//
// The metric is stationary and axisymmetric, so it depends only on r and th;
// d/dt and d/dphi of every component vanish identically.  That is exactly why
// E = -p_t and L_z = p_phi are constants of the motion, and it lets the
// geodesic integrator carry those two momenta exactly.
// ---------------------------------------------------------------------------
struct MetricUpper {
    double tt, tp, rr, thth, pp;
};

inline MetricUpper metric_upper(const Geom& g) {
    MetricUpper m;
    const double SD = g.Sigma * g.Delta;
    m.tt   = -g.A / SD;
    m.tp   = -2.0 * g.a * g.r / SD;
    m.rr   = g.Delta / g.Sigma;
    m.thth = 1.0 / g.Sigma;
    m.pp   = (g.Delta - g.a2 * g.sn2) / (SD * g.sn2);
    return m;
}

inline MetricUpper metric_upper(double a, double r, double th) {
    return metric_upper(geom_at(a, r, th));
}

// Derivatives of g^{mu nu} with respect to r (d_r) and theta (d_th).
struct MetricUpperDeriv {
    MetricUpper g;      // the metric itself
    MetricUpper dr;     // d/dr
    MetricUpper dth;    // d/dtheta
};

inline MetricUpperDeriv metric_upper_deriv(const Geom& g) {
    MetricUpperDeriv o;
    o.g = metric_upper(g);

    const double r = g.r, a2 = g.a2, a = g.a;
    const double sn = g.sn, cs = g.cs, sn2 = g.sn2;

    // Building blocks and their derivatives.
    const double S = g.Sigma, D = g.Delta, A = g.A;
    const double S_r = 2.0 * r,            S_th = -2.0 * a2 * sn * cs;
    const double D_r = 2.0 * r - 2.0,      D_th = 0.0;
    const double A_r = 4.0 * r * (r * r + a2) - a2 * D_r * sn2;
    const double A_th = -2.0 * a2 * D * sn * cs;

    const double Q   = S * D;                       // common denominator
    const double Q_r = S_r * D + S * D_r;
    const double Q_th = S_th * D + S * D_th;
    const double Q2  = Q * Q;

    // g^tt = -A / (S D)
    o.dr.tt  = -(A_r * Q - A * Q_r) / Q2;
    o.dth.tt = -(A_th * Q - A * Q_th) / Q2;

    // g^tphi = -2 a r / (S D)
    const double N = -2.0 * a * r, N_r = -2.0 * a, N_th = 0.0;
    o.dr.tp  = (N_r * Q - N * Q_r) / Q2;
    o.dth.tp = (N_th * Q - N * Q_th) / Q2;

    // g^rr = D / S
    o.dr.rr  = (D_r * S - D * S_r) / (S * S);
    o.dth.rr = (D_th * S - D * S_th) / (S * S);

    // g^thth = 1 / S
    o.dr.thth  = -S_r / (S * S);
    o.dth.thth = -S_th / (S * S);

    // g^phiphi = (D - a^2 sin^2) / (S D sin^2)
    const double P    = D - a2 * sn2;
    const double P_r  = D_r;
    const double P_th = -2.0 * a2 * sn * cs;
    const double R    = Q * sn2;
    const double R_r  = Q_r * sn2;
    const double R_th = Q_th * sn2 + Q * 2.0 * sn * cs;
    o.dr.pp  = (P_r * R - P * R_r) / (R * R);
    o.dth.pp = (P_th * R - P * R_th) / (R * R);

    return o;
}

inline MetricUpperDeriv metric_upper_deriv(double a, double r, double th) {
    return metric_upper_deriv(geom_at(a, r, th));
}

// ---------------------------------------------------------------------------
// Horizons, ergosphere, and the special circular orbits.
// ---------------------------------------------------------------------------

// Outer event horizon r_+ = M + sqrt(M^2 - a^2).
inline double horizon_outer(double a) {
    const double d = 1.0 - a * a;
    return 1.0 + std::sqrt(std::max(0.0, d));
}
inline double horizon_inner(double a) {
    const double d = 1.0 - a * a;
    return 1.0 - std::sqrt(std::max(0.0, d));
}

// Angular velocity of the horizon itself, Omega_H = a c^3 / (2 G M r_+).
//
// A Kerr horizon is not a surface that can be "at rest": it drags spacetime
// round with it rigidly, at this rate.  Every zero-angular-momentum observer
// approaches Omega_H as they approach the horizon, which is checked in
// validate.cpp.  For a hole at the Thorne limit it is 0.469 c^3/GM - so a
// billion-solar-mass horizon, sixty-six astronomical units around, turns once
// in eighteen hours.
inline double horizon_angular_velocity(double a) {
    const double rp = horizon_outer(a);
    return (rp > 0.0) ? a / (2.0 * rp) : 0.0;
}

// Static limit (outer boundary of the ergosphere): r_E = M + sqrt(M^2-a^2cos^2).
inline double ergosphere_outer(double a, double th) {
    const double c = std::cos(th);
    return 1.0 + std::sqrt(std::max(0.0, 1.0 - a * a * c * c));
}

// Equatorial circular photon orbit (Bardeen, Press & Teukolsky 1972).
// sigma = +1 prograde, -1 retrograde.
inline double photon_circular_orbit(double a, int sigma) {
    const double x = std::clamp(-static_cast<double>(sigma) * a, -1.0, 1.0);
    return 2.0 * (1.0 + std::cos((2.0 / 3.0) * std::acos(x)));
}

// Innermost stable circular orbit, Bardeen-Press-Teukolsky (1972) eq. 2.21.
inline double isco_radius(double a, int sigma = +1) {
    const double a2 = a * a;
    const double Z1 = 1.0 + std::cbrt(1.0 - a2) * (std::cbrt(1.0 + a) + std::cbrt(1.0 - a));
    const double Z2 = std::sqrt(3.0 * a2 + Z1 * Z1);
    const double s  = std::sqrt(std::max(0.0, (3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2)));
    return 3.0 + Z2 - static_cast<double>(sigma) * s;
}

// Specific energy E = -u_t of an equatorial circular geodesic at radius r.
inline double circular_energy(double a, double r, int sigma = +1) {
    const double sr = std::sqrt(r);
    const double num = r * r - 2.0 * r + sigma * a * sr;
    const double den = r * std::sqrt(r * r - 3.0 * r + 2.0 * sigma * a * sr);
    return num / den;
}

// Specific angular momentum L_z = u_phi of an equatorial circular geodesic.
inline double circular_angmom(double a, double r, int sigma = +1) {
    const double sr = std::sqrt(r);
    const double num = sigma * sr * (r * r - 2.0 * sigma * a * sr + a * a);
    const double den = r * std::sqrt(r * r - 3.0 * r + 2.0 * sigma * a * sr);
    return num / den;
}

// Radiative efficiency of a Novikov-Thorne disc: the binding energy per unit
// rest mass released between infinity and the ISCO.  0.0572 for a = 0,
// 0.42 for a -> 1.
inline double disc_efficiency(double a, int sigma = +1) {
    return 1.0 - circular_energy(a, isco_radius(a, sigma), sigma);
}

// ---------------------------------------------------------------------------
// Observers and orbiting matter.
// ---------------------------------------------------------------------------

// Four-velocity (u^t, u^r, u^th, u^phi) of a circular Keplerian orbit in the
// equatorial plane.  This is the motion of the disc material.
inline std::array<double, 4> keplerian_four_velocity(double a, double r, int sigma = +1) {
    const double r32 = r * std::sqrt(r);
    const double Omega = sigma / (r32 + sigma * a);
    const double disc  = r * r * r - 3.0 * r * r + 2.0 * sigma * a * r32;
    const double ut    = (r32 + sigma * a) / std::sqrt(std::max(1e-300, disc));
    return {ut, 0.0, 0.0, Omega * ut};
}

// Angular velocity of a zero-angular-momentum observer (ZAMO):
//   omega = -g_{t phi} / g_{phi phi} = 2 M a r / A.
inline double zamo_omega(const Geom& g) { return 2.0 * g.a * g.r / g.A; }

// Lapse of the ZAMO frame: alpha = sqrt(Sigma Delta / A).
inline double zamo_lapse(const Geom& g) { return std::sqrt(g.Sigma * g.Delta / g.A); }

// Orthonormal tetrad e_(a)^mu of the locally non-rotating (ZAMO) frame.
// Rows: 0 = timelike, 1 = radial, 2 = polar, 3 = azimuthal.
// This frame is regular everywhere outside the horizon, including inside the
// ergosphere where no static observer exists.
struct Tetrad {
    double e[4][4];  // e[a][mu]
};

inline Tetrad zamo_tetrad(const Geom& g) {
    Tetrad T{};
    const double alpha = zamo_lapse(g);
    const double omega = zamo_omega(g);
    T.e[0][0] = 1.0 / alpha;
    T.e[0][3] = omega / alpha;
    T.e[1][1] = std::sqrt(g.Delta / g.Sigma);
    T.e[2][2] = 1.0 / std::sqrt(g.Sigma);
    T.e[3][3] = std::sqrt(g.Sigma / g.A) / g.sn;
    return T;
}

// ---------------------------------------------------------------------------
// Contractions used everywhere.
// ---------------------------------------------------------------------------

// Hamiltonian contraction  H = g^{mu nu} p_mu p_nu.  Zero for photons,
// -m^2 (= -1 for unit rest mass) for massive particles.
inline double hamiltonian(const MetricUpper& m, const std::array<double, 4>& p) {
    return m.tt * p[0] * p[0]
         + 2.0 * m.tp * p[0] * p[3]
         + m.rr * p[1] * p[1]
         + m.thth * p[2] * p[2]
         + m.pp * p[3] * p[3];
}

// Sum of the magnitudes of the individual terms in that contraction.  Dividing
// by it turns the null condition into a scale-free measure of how much
// cancellation actually survived, which is the meaningful accuracy statement:
// close to the horizon the individual terms grow like 1/Delta, so a fixed
// absolute threshold on H would be far too strict there and far too lax far
// away.
inline double hamiltonian_scale(const MetricUpper& m, const std::array<double, 4>& p) {
    return std::fabs(m.tt * p[0] * p[0])
         + std::fabs(2.0 * m.tp * p[0] * p[3])
         + std::fabs(m.rr * p[1] * p[1])
         + std::fabs(m.thth * p[2] * p[2])
         + std::fabs(m.pp * p[3] * p[3]);
}

// Raise an index: p^mu = g^{mu nu} p_nu.
inline std::array<double, 4> raise(const MetricUpper& m, const std::array<double, 4>& p) {
    return {m.tt * p[0] + m.tp * p[3],
            m.rr * p[1],
            m.thth * p[2],
            m.tp * p[0] + m.pp * p[3]};
}

// Lower an index: p_mu = g_{mu nu} p^nu.
inline std::array<double, 4> lower(const MetricLower& m, const std::array<double, 4>& u) {
    return {m.tt * u[0] + m.tp * u[3],
            m.rr * u[1],
            m.thth * u[2],
            m.tp * u[0] + m.pp * u[3]};
}

// Carter's constant, the "fourth" integral that makes Kerr geodesics
// separable (Carter 1968):
//   Q = p_th^2 + cos^2(th) [ a^2 (m^2 - E^2) + L_z^2 / sin^2(th) ].
inline double carter_constant(double a, double th, double p_th, double E, double Lz,
                              double m2 = 0.0) {
    double sn = std::sin(th);
    if (std::fabs(sn) < kSinThetaMin) sn = kSinThetaMin;
    const double cs2 = std::cos(th) * std::cos(th);
    return p_th * p_th + cs2 * (a * a * (m2 - E * E) + Lz * Lz / (sn * sn));
}

} // namespace bh
