// tov.h - Hydrostatic equilibrium in general relativity.
//
// The Tolman-Oppenheimer-Volkoff equation is what Einstein's field equations
// reduce to for a static, spherically symmetric star.  Writing the interior
// metric as
//
//     ds^2 = -e^{2 Phi(r)} c^2 dt^2 + (1 - 2 G m(r) / (r c^2))^{-1} dr^2
//            + r^2 dOmega^2
//
// the (tt) and (rr) components of G_{mu nu} = 8 pi T_{mu nu} / c^4 for a
// perfect fluid give
//
//     dm/dr = 4 pi r^2 rho
//
//     dP/dr = - G (rho + P/c^2)(m + 4 pi r^3 P / c^2)
//             ---------------------------------------
//                    r^2 (1 - 2 G m / (r c^2))
//
//     dPhi/dr = -1/(rho c^2 + P) dP/dr
//
// Every difference from Newtonian hydrostatics is visible in that second
// equation: pressure gravitates (the P/c^2 terms), and the denominator blows
// up as the star approaches its own Schwarzschild radius.  Both push the same
// way - towards collapse - which is why a relativistic star has a *maximum*
// mass at all, and why no equation of state can evade one.
//
// Integrating outward from a chosen central density until P reaches zero gives
// one star: its radius R and its gravitational mass M = m(R).  Sweeping the
// central density traces out the mass-radius curve of the equation of state,
// whose turning point is the maximum mass.
#pragma once

#include "eos.h"
#include "constants.h"

#include <cmath>
#include <vector>
#include <algorithm>

namespace ns {

struct Star {
    double R = 0.0;              // circumferential radius, m
    double M = 0.0;              // gravitational mass, kg
    double rho_c = 0.0;          // central density, kg m^-3
    double P_c = 0.0;            // central pressure, Pa
    double compactness = 0.0;    // r_s / R = 2GM/(Rc^2)
    double redshift = 0.0;       // surface gravitational redshift z
    double max_sound_speed = 0.0;// max c_s/c inside; > 1 means acausal
    double binding_energy = 0.0; // (baryon mass - gravitational mass) c^2, J
    double baryon_mass = 0.0;    // kg
    bool   ok = false;

    double M_solar() const { return M / phys::M_sun; }
    double R_km() const { return R / 1000.0; }
};

// Integrate one star from a given central density.
//
// The scheme is fourth-order Runge-Kutta in r with an adaptive step tied to
// the pressure scale height, which matters because P falls by many orders of
// magnitude over the last few hundred metres of the crust.  The surface is
// located by bisection on P = P_surface rather than by whichever step happens
// to overshoot.
inline Star solve_tov(const EOS& eos, double rho_c,
                      double P_surface_frac = 1e-10,
                      int max_steps = 2000000) {
    using phys::G;
    using phys::c;
    const double c2 = c * c;

    Star s;
    s.rho_c = rho_c;
    s.P_c = eos.pressure(rho_c);
    if (!(s.P_c > 0.0)) return s;

    const double P_stop = s.P_c * P_surface_frac;

    // State: y = (m, P, Phi, baryon-mass integral)
    // Phi is fixed up at the end by matching to the exterior Schwarzschild
    // solution, so its starting value is arbitrary.
    struct Y { double m, P, phi, mb; };

    auto deriv = [&](double r, const Y& y) -> Y {
        Y d{};
        if (r <= 0.0) return d;
        const double rho = eos.density(y.P);
        const double f = 1.0 - 2.0 * G * y.m / (r * c2);
        d.m = 4.0 * M_PI * r * r * rho;
        if (f <= 1e-12) {         // inside its own horizon: not a star
            d.P = -1e300;
            d.phi = 0.0;
            d.mb = 0.0;
            return d;
        }
        d.P = -G * (rho + y.P / c2) * (y.m + 4.0 * M_PI * r * r * r * y.P / c2) /
              (r * r * f);
        d.phi = -d.P / (rho * c2 + y.P);
        // Proper volume element for the baryon (rest) mass.
        d.mb = 4.0 * M_PI * r * r * rho / std::sqrt(f);
        return d;
    };

    // Series expansion near the centre avoids the 0/0 at r = 0:
    //   m  ~ (4/3) pi rho_c r^3
    //   P  ~ P_c - (2 pi G / 3)(rho_c + P_c/c^2)(rho_c + 3 P_c/c^2) r^2
    double r = 1e-3;
    Y y{};
    y.m = (4.0 / 3.0) * M_PI * rho_c * r * r * r;
    y.P = s.P_c - (2.0 * M_PI * G / 3.0) * (rho_c + s.P_c / c2) *
                  (rho_c + 3.0 * s.P_c / c2) * r * r;
    y.phi = 0.0;
    y.mb = 0.0;

    s.max_sound_speed = std::sqrt(std::max(0.0, eos.sound_speed_squared(rho_c)));

    double dr = 1.0;   // metres
    for (int step = 0; step < max_steps; ++step) {
        // Step size from the local pressure scale height, clamped so the
        // integration neither crawls nor jumps over the surface.
        const Y d0 = deriv(r, y);
        if (d0.P < -1e299) return s;
        const double scale = (d0.P != 0.0) ? std::fabs(y.P / d0.P) : 1e30;
        dr = std::clamp(0.02 * scale, 0.1, 200.0);

        auto add = [](const Y& a, const Y& b, double h) {
            return Y{a.m + h * b.m, a.P + h * b.P, a.phi + h * b.phi, a.mb + h * b.mb};
        };
        const Y k1 = d0;
        const Y k2 = deriv(r + 0.5 * dr, add(y, k1, 0.5 * dr));
        const Y k3 = deriv(r + 0.5 * dr, add(y, k2, 0.5 * dr));
        const Y k4 = deriv(r + dr, add(y, k3, dr));

        Y yn{};
        yn.m   = y.m   + dr / 6.0 * (k1.m + 2 * k2.m + 2 * k3.m + k4.m);
        yn.P   = y.P   + dr / 6.0 * (k1.P + 2 * k2.P + 2 * k3.P + k4.P);
        yn.phi = y.phi + dr / 6.0 * (k1.phi + 2 * k2.phi + 2 * k3.phi + k4.phi);
        yn.mb  = y.mb  + dr / 6.0 * (k1.mb + 2 * k2.mb + 2 * k3.mb + k4.mb);

        s.max_sound_speed = std::max(
            s.max_sound_speed,
            std::sqrt(std::max(0.0, eos.sound_speed_squared(eos.density(y.P)))));

        if (yn.P <= P_stop || !std::isfinite(yn.P)) {
            // Bisect on r within this step to land on the surface.
            double lo = 0.0, hi = dr;
            for (int k = 0; k < 80; ++k) {
                const double mid = 0.5 * (lo + hi);
                const Y a1 = deriv(r, y);
                const Y a2 = deriv(r + 0.5 * mid, add(y, a1, 0.5 * mid));
                const Y a3 = deriv(r + 0.5 * mid, add(y, a2, 0.5 * mid));
                const Y a4 = deriv(r + mid, add(y, a3, mid));
                const double Pm = y.P + mid / 6.0 * (a1.P + 2 * a2.P + 2 * a3.P + a4.P);
                if (Pm <= P_stop || !std::isfinite(Pm)) hi = mid; else lo = mid;
            }
            const double h = 0.5 * (lo + hi);
            const Y a1 = deriv(r, y);
            const Y a2 = deriv(r + 0.5 * h, add(y, a1, 0.5 * h));
            const Y a3 = deriv(r + 0.5 * h, add(y, a2, 0.5 * h));
            const Y a4 = deriv(r + h, add(y, a3, h));
            s.R = r + h;
            s.M = y.m + h / 6.0 * (a1.m + 2 * a2.m + 2 * a3.m + a4.m);
            s.baryon_mass = y.mb + h / 6.0 * (a1.mb + 2 * a2.mb + 2 * a3.mb + a4.mb);
            s.ok = s.R > 0.0 && s.M > 0.0;
            break;
        }
        y = yn;
        r += dr;
    }

    if (!s.ok) return s;
    s.compactness = 2.0 * G * s.M / (s.R * c2);
    if (s.compactness < 1.0)
        s.redshift = 1.0 / std::sqrt(1.0 - s.compactness) - 1.0;
    s.binding_energy = (s.baryon_mass - s.M) * c2;
    return s;
}

// Sweep central density to trace the mass-radius curve.  The maximum of M
// along it is the maximum mass the equation of state can support; stars beyond
// that turning point are unstable to radial oscillation and collapse.
struct MassRadiusCurve {
    std::vector<Star> stars;
    Star max_mass;
    double M_max_solar() const { return max_mass.M / phys::M_sun; }
};

inline MassRadiusCurve mass_radius_curve(const EOS& eos, double rho_lo, double rho_hi,
                                         int n = 60) {
    MassRadiusCurve curve;
    curve.stars.reserve(n);
    for (int i = 0; i < n; ++i) {
        const double f = static_cast<double>(i) / (n - 1);
        const double rho = rho_lo * std::pow(rho_hi / rho_lo, f);
        const Star s = solve_tov(eos, rho);
        if (!s.ok) continue;
        curve.stars.push_back(s);
        if (s.M > curve.max_mass.M) curve.max_mass = s;
    }
    return curve;
}

// Find the central density that produces a star of a given gravitational mass,
// on the stable branch (below the maximum-mass turning point).
inline Star star_of_mass(const EOS& eos, double M_target,
                         double rho_lo = 2e17, double rho_hi = 5e18) {
    // M(rho_c) increases monotonically on the stable branch, so bisect there.
    const MassRadiusCurve curve = mass_radius_curve(eos, rho_lo, rho_hi, 40);
    if (curve.max_mass.M < M_target) return Star{};   // this EOS cannot do it

    double lo = rho_lo, hi = curve.max_mass.rho_c;
    for (int i = 0; i < 80; ++i) {
        const double mid = std::sqrt(lo * hi);
        const Star s = solve_tov(eos, mid);
        if (!s.ok || s.M < M_target) lo = mid; else hi = mid;
    }
    return solve_tov(eos, std::sqrt(lo * hi));
}

// ---------------------------------------------------------------------------
// The exact interior Schwarzschild solution, for validation.
//
// For rho = const the TOV equation integrates in closed form (Schwarzschild
// 1916).  With y(r) = sqrt(1 - 2GM r^2/(R^3 c^2)) and y_R = sqrt(1 - r_s/R),
//
//     P(r) = rho c^2 (y(r) - y_R) / (3 y_R - y(r))
//
// and the central pressure diverges when y_R = 1/3, i.e. when the compactness
// reaches 8/9.  That is Buchdahl's bound: no static star of any equation of
// state can be more compact than R = 9 GM / (4 c^2).
// ---------------------------------------------------------------------------
inline double interior_schwarzschild_pressure(double rho, double R, double r) {
    const double M = (4.0 / 3.0) * M_PI * R * R * R * rho;
    const double rs = 2.0 * phys::G * M / (phys::c * phys::c);
    const double yR = std::sqrt(std::max(0.0, 1.0 - rs / R));
    const double y  = std::sqrt(std::max(0.0, 1.0 - rs * r * r / (R * R * R)));
    const double den = 3.0 * yR - y;
    if (std::fabs(den) < 1e-30) return 1e300;
    return rho * phys::c * phys::c * (y - yR) / den;
}

// Buchdahl's bound on compactness, 8/9.
inline constexpr double kBuchdahlCompactness = 8.0 / 9.0;

} // namespace ns
