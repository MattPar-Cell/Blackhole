// neutronstar.h - Rendering a neutron star surface.
//
// A neutron star is a *body*, not a hole, and that changes what the picture
// shows.  There is no horizon and no shadow; light stops at a surface.  For a
// 1.4 solar mass star of radius 12 km the surface sits at r = 5.8 GM/c^2,
// comfortably outside the photon sphere at 3 GM/c^2, so no light can orbit and
// there is no photon ring either.
//
// What there is instead is extreme lensing of the star's own surface.  Rays
// that leave tangentially still curve round to the observer, so a large slice
// of the far hemisphere is visible: for typical compactness about 76% of the
// total surface area, not 50%.  Beloborodov (2002) gives the accurate
// approximation
//
//     1 - cos(alpha) = (1 - cos(psi)) (1 - r_s/R)
//
// relating the emission angle alpha from the local normal to the angle psi
// between the emission point and the line of sight; setting alpha = 90 degrees
// gives the last visible point.  validate.cpp checks the traced geodesics
// against both that approximation and an exact quadrature.
//
// Spacetime: by Birkhoff's theorem the exterior of a static spherical star is
// *exactly* Schwarzschild, so the same metric and integrator used for the black
// holes apply unchanged, with a = 0.  Rotation is handled the way NICER's
// pulse-profile modelling does in its "Schwarzschild + Doppler" form: the
// metric stays Schwarzschild while the surface moves, giving the Doppler shift,
// aberration and beaming of a rotating surface.  What that leaves out is the
// frame dragging and the mass quadrupole of a genuinely rotating star - Kerr is
// *not* the exterior of a rotating neutron star, because its quadrupole moment
// is not the one a real star has.  The approximation is good to a few percent
// below a few hundred hertz and degrades for the fastest millisecond pulsars.
#pragma once

#include "kerr.h"
#include "constants.h"
#include "tov.h"

#include <cmath>
#include <array>

namespace bh {

struct NeutronStar {
    bool   enabled = false;

    // Structure, from a TOV solve.
    double M_kg = 1.4 * phys::M_sun;
    double R_m  = 12000.0;
    double R    = 5.8;          // surface radius in GM/c^2

    // Surface emission.
    double T_eff = 1.0e6;       // K, effective temperature of the bulk surface

    // Rotation.  spin_hz is the observed spin frequency; Omega is in geometric
    // units (c^3/GM), which is what the metric wants.
    double spin_hz = 0.0;
    double Omega   = 0.0;

    // Magnetic polar caps.  The magnetic axis is tilted by `cap_tilt` from the
    // spin axis and carried round with the star, which is what makes a pulsar
    // pulse.  Two antipodal caps, as for a dipole field.
    bool   caps = false;
    double cap_tilt = 60.0 * M_PI / 180.0;   // alpha, magnetic obliquity
    double cap_radius = 20.0 * M_PI / 180.0; // angular radius of a cap
    double cap_T = 3.0e6;                    // K
    double cap_phase0 = 0.0;
    bool   two_caps = true;

    void set_from(const ns::Star& s) {
        M_kg = s.M;
        R_m  = s.R;
        R    = s.R / phys::r_g_metres(s.M);
    }

    void set_spin(double hz) {
        spin_hz = hz;
        // Omega in geometric units: multiply by the light-crossing time GM/c^3.
        Omega = 2.0 * M_PI * hz * phys::r_g_metres(M_kg) / phys::c;
    }

    // Equatorial surface speed as measured by a local static observer, in c.
    double equatorial_speed() const {
        const double f = 1.0 - 2.0 / R;
        return (f > 0.0) ? Omega * R / std::sqrt(f) : 0.0;
    }

    // Four-velocity of the rotating surface at colatitude theta.  Exact for
    // rigid rotation in Schwarzschild: the only non-zero components are t and
    // phi, normalised by u.u = -1.
    std::array<double, 4> surface_four_velocity(double theta) const {
        const double sn = std::sin(theta);
        const double norm2 = 1.0 - 2.0 / R - Omega * Omega * R * R * sn * sn;
        const double ut = (norm2 > 1e-12) ? 1.0 / std::sqrt(norm2) : 1.0;
        return {ut, 0.0, 0.0, Omega * ut};
    }

    // Local surface temperature at a point, allowing for the polar caps.
    // `t_emit` is the coordinate time the light left, so the caps are found
    // where they were then rather than where they are now.
    double temperature_at(double theta, double phi, double t_emit) const {
        if (!caps) return T_eff;

        const double spin_phase = cap_phase0 + Omega * t_emit;
        // Magnetic axis, tilted from the spin axis and carried round with it.
        const double m[3] = {std::sin(cap_tilt) * std::cos(spin_phase),
                             std::sin(cap_tilt) * std::sin(spin_phase),
                             std::cos(cap_tilt)};
        const double n[3] = {std::sin(theta) * std::cos(phi),
                             std::sin(theta) * std::sin(phi),
                             std::cos(theta)};
        const double dot = n[0] * m[0] + n[1] * m[1] + n[2] * m[2];
        double ang = std::acos(std::clamp(dot, -1.0, 1.0));
        if (two_caps) ang = std::min(ang, M_PI - ang);   // the antipodal cap

        if (ang <= cap_radius) {
            // Soften the rim over the outer tenth so it does not alias; the
            // real boundary of an open field line region is not a step either.
            const double edge = 0.9 * cap_radius;
            if (ang <= edge) return cap_T;
            const double u = (ang - edge) / (cap_radius - edge);
            const double w = 0.5 * (1.0 + std::cos(M_PI * u));   // smooth to 0
            return T_eff + (cap_T - T_eff) * w;
        }
        return T_eff;
    }

    // Bolometric luminosity of the bulk surface as seen from far away, in W.
    // The proper area is 4 pi R^2 (the circumferential radius is defined that
    // way) and the redshift takes two factors of (1 - r_s/R).
    double luminosity_infinity() const {
        const double f = 1.0 - 2.0 / R;
        return 4.0 * M_PI * R_m * R_m * phys::sigma_SB *
               std::pow(T_eff, 4.0) * f;
    }

    double surface_redshift() const {
        const double f = 1.0 - 2.0 / R;
        return (f > 0.0) ? 1.0 / std::sqrt(f) - 1.0 : 0.0;
    }

    // Fraction of the total surface an observer at infinity can see, using
    // Beloborodov's relation with alpha = 90 degrees.
    double visible_fraction() const {
        const double rs_over_R = 2.0 / R;
        const double cos_psi = 1.0 - 1.0 / (1.0 - rs_over_R);
        return 0.5 * (1.0 - std::clamp(cos_psi, -1.0, 1.0));
    }
};

// Exact maximum visible colatitude, by quadrature rather than approximation.
//
// A photon leaving the surface tangentially has impact parameter
// b = R / sqrt(1 - r_s/R); the angle it sweeps on its way out to infinity is
//
//     psi = Integral_R^inf  (b / r^2) dr / sqrt(1 - b^2 (1 - r_s/r) / r^2)
//
// The integrand has an inverse-square-root singularity at r = R, removed here
// by substituting r = R / (1 - s^2).
inline double max_visible_angle_exact(double R_geo) {
    const double rs = 2.0;
    const double b = R_geo / std::sqrt(1.0 - rs / R_geo);

    auto integrand = [&](double s) {
        // r = R/(1-s^2), dr = 2 R s /(1-s^2)^2 ds
        const double one_minus = 1.0 - s * s;
        if (one_minus <= 1e-14) return 0.0;
        const double r = R_geo / one_minus;
        const double dr_ds = 2.0 * R_geo * s / (one_minus * one_minus);
        const double f = 1.0 - b * b * (1.0 - rs / r) / (r * r);
        if (f <= 0.0) return 0.0;
        return (b / (r * r)) * dr_ds / std::sqrt(f);
    };

    // Composite Gauss-Legendre over s in [0, 1).
    static const double gx[8] = {-0.9602898564975363, -0.7966664774136267,
                                 -0.5255324099163290, -0.1834346424956498,
                                  0.1834346424956498,  0.5255324099163290,
                                  0.7966664774136267,  0.9602898564975363};
    static const double gw[8] = {0.1012285362903763, 0.2223810344533745,
                                 0.3137066458778873, 0.3626837833783620,
                                 0.3626837833783620, 0.3137066458778873,
                                 0.2223810344533745, 0.1012285362903763};
    const int panels = 20000;
    const double s_max = 1.0 - 1e-9;
    double total = 0.0;
    for (int p = 0; p < panels; ++p) {
        const double a0 = s_max * p / panels, a1 = s_max * (p + 1) / panels;
        const double hm = 0.5 * (a1 - a0), cm = 0.5 * (a1 + a0);
        for (int k = 0; k < 8; ++k) total += hm * gw[k] * integrand(cm + hm * gx[k]);
    }
    return total;
}

} // namespace bh
