// disc.h - Relativistic thin accretion disc (Novikov & Thorne 1973,
//          Page & Thorne 1974).
//
// The disc is geometrically thin, optically thick, and in the equatorial
// plane.  Each annulus orbits on a circular Keplerian geodesic, radiates as a
// blackbody from its own rest frame, and the inner edge is the ISCO with a
// zero-torque boundary condition.
//
// The radiated flux comes from conservation of rest mass, angular momentum and
// energy in the disc, which give (Page & Thorne 1974, eq. 15):
//
//                 Mdot        -dOmega/dr      r
//     F(r)  =  ---------- . -------------- . I  (E - Omega L) dL/dr dr'
//               4 pi r        (E - Omega L)^2  r_isco
//
// with E(r), L(r), Omega(r) the specific energy, specific angular momentum and
// angular velocity of a circular equatorial Kerr geodesic.  Page & Thorne
// evaluated that integral in closed form; the result is `page_thorne_W` below.
//
// validate.cpp re-derives the same quantity by brute-force numerical
// quadrature of the integral above and checks the two agree (they do, to about
// 1e-8 over four decades in radius and four values of the spin), then checks
// the two things the result must independently satisfy:
//
//   * the far-field limit  F -> 3 G M Mdot / (8 pi r^3), the Shakura-Sunyaev
//     law, and
//   * the energy budget  Integral F E dA = eta Mdot c^2, eta = 1 - E(r_isco),
//     where the factor E(r) redshifts locally emitted energy to infinity.
//
// Effective temperature follows from F = sigma T^4 (one face).
#pragma once

#include "kerr.h"
#include "constants.h"

#include <cmath>
#include <algorithm>

namespace bh {

// Page & Thorne's closed-form evaluation of
//     W(r) = Integral_{r_isco}^{r} (E - Omega L) dL/dr' dr'
// in terms of x = sqrt(r/M) and the three roots of x^3 - 3x + 2a = 0, which
// are the square roots of the radii where the circular-orbit denominators
// vanish (x1 is the prograde photon orbit).
inline double page_thorne_W(double a, double r, double r_isco) {
    const double x  = std::sqrt(r);
    const double x0 = std::sqrt(r_isco);
    const double ac = std::acos(std::clamp(a, -1.0, 1.0));
    const double x1 =  2.0 * std::cos(ac / 3.0 - M_PI / 3.0);
    const double x2 =  2.0 * std::cos(ac / 3.0 + M_PI / 3.0);
    const double x3 = -2.0 * std::cos(ac / 3.0);

    auto term = [&](double xi, double xj, double xk) {
        // As a -> 0 the middle root goes to zero and so does its coefficient
        // (numerator (x_i - a)^2 -> 0 faster than x_i); the limit is exactly 0.
        if (std::fabs(xi) < 1e-13) return 0.0;
        return 3.0 * (xi - a) * (xi - a) / (xi * (xi - xj) * (xi - xk)) *
               std::log((x - xi) / (x0 - xi));
    };

    return x - x0 - 1.5 * a * std::log(x / x0)
         - term(x1, x2, x3) - term(x2, x1, x3) - term(x3, x1, x2);
}

class NovikovThorneDisc {
public:
    // M_kg      : black hole mass
    // spin      : dimensionless a
    // r_out     : outer radius, in r_g
    // eddington : accretion rate as a fraction of the Eddington rate
    // prograde  : disc co-rotates (+1) or counter-rotates (-1) with the hole
    NovikovThorneDisc(double M_kg, double spin, double r_out, double eddington,
                      int prograde = +1)
        : M_(M_kg), a_(std::clamp(spin, -0.9999, 0.9999)), sigma_(prograde),
          r_in_(isco_radius(std::clamp(spin, -0.9999, 0.9999), prograde)),
          r_out_(std::max(r_out, isco_radius(std::clamp(spin, -0.9999, 0.9999), prograde) * 1.02)) {
        eta_   = disc_efficiency(a_, sigma_);
        mdot_  = eddington * phys::mdot_eddington(M_kg, eta_);   // kg/s
        r_g_   = phys::r_g_metres(M_kg);
        build_table();
    }

    double inner_radius() const { return r_in_; }
    double outer_radius() const { return r_out_; }
    double efficiency()   const { return eta_; }
    double mdot_si()      const { return mdot_; }
    double luminosity()   const { return eta_ * mdot_ * phys::c * phys::c; }  // W
    int    sense()        const { return sigma_; }

    bool contains(double r) const { return r >= r_in_ && r <= r_out_; }

    // Flux from one face, W m^-2, as measured in the local rest frame.
    double flux(double r) const {
        if (r <= r_in_ || r > r_out_) return 0.0;
        return flux_dimensionless(r) * mdot_ * phys::c * phys::c / (r_g_ * r_g_);
    }

    // Effective blackbody temperature in the local rest frame, in kelvin.
    double temperature(double r) const {
        const double F = flux(r);
        if (F <= 0.0) return 0.0;
        return std::pow(F / phys::sigma_SB, 0.25);
    }

    // Peak temperature over the disc, useful for exposure setting.
    double peak_temperature() const { return T_peak_; }

    // Dimensionless flux F/(Mdot c^2 / r_g^2) with r in gravitational radii.
    double flux_dimensionless(double r) const {
        if (r <= r_in_) return 0.0;
        const double W = page_thorne_W(a_, r, r_in_);
        const double Om  = omega(r);
        const double dOm = omega_prime(r);
        const double E   = circular_energy(a_, r, sigma_);
        const double L   = circular_angmom(a_, r, sigma_);
        const double den = (E - Om * L);
        if (!(den > 0.0)) return 0.0;
        const double F = (-dOm) / (4.0 * M_PI * r * den * den) * W;
        return (F > 0.0 && std::isfinite(F)) ? F : 0.0;
    }

    // Four-velocity of the disc material at radius r (equatorial).
    std::array<double, 4> four_velocity(double r) const {
        return keplerian_four_velocity(a_, r, sigma_);
    }

private:
    double M_, a_;
    int    sigma_;
    double r_in_, r_out_;
    double eta_ = 0.0, mdot_ = 0.0, r_g_ = 1.0, T_peak_ = 0.0;

    double omega(double r) const {
        const double r32 = r * std::sqrt(r);
        return sigma_ / (r32 + sigma_ * a_);
    }
    double omega_prime(double r) const {
        // d/dr [ sigma / (r^{3/2} + sigma a) ] = -sigma (3/2) sqrt(r) / (...)^2
        const double r32 = r * std::sqrt(r);
        const double d   = r32 + sigma_ * a_;
        return -sigma_ * 1.5 * std::sqrt(r) / (d * d);
    }
    void build_table() {
        // Peak temperature, used only to pick sensible default exposures.
        for (int i = 0; i < 4096; ++i) {
            const double r = r_in_ * std::pow(r_out_ / r_in_, i / 4095.0);
            T_peak_ = std::max(T_peak_, temperature(r));
        }
    }
};

} // namespace bh
