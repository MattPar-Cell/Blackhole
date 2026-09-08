// jet.h - Relativistic plasma jets along the spin axis.
//
// A quasar's jet is not powered by the accretion disc.  It is powered by the
// *rotation of the black hole itself*.  Magnetic field lines threading the
// horizon are dragged round with it, and because a rotating horizon behaves
// like a resistive membrane the twist propagates outwards as a Poynting flux.
// Blandford & Znajek (1977) derived the power; the modern calibrated form,
// from general-relativistic magnetohydrodynamic simulations of the
// magnetically arrested state (Tchekhovskoy, Narayan & McKinney 2011), is
//
//     P  =  kappa/(4 pi) . Phi_H^2 . Omega_H^2 . f(Omega_H),
//     f  =  1 + 1.38 Omega_H^2 - 9.2 Omega_H^4,        kappa = 0.053
//
// with Phi_H = phi sqrt(Mdot r_g^2 c) the magnetic flux threading the horizon
// and phi ~ 50 the saturation value of a magnetically arrested disc.  Written
// as a fraction of the accreted rest-mass energy this is
//
//     P / (Mdot c^2)  =  kappa/(4 pi) . phi^2 . Omega_H^2 . f(Omega_H)
//
// which for the Thorne-limit hole of the `quasar` preset (a = 0.998,
// Omega_H = 0.4693) comes to 1.99.  The jet carries away *twice* the energy
// that the accreted matter brought in.  That is not a bookkeeping error: the
// surplus is extracted from the hole's rotational energy, so a jet at this
// spin is actively spinning the hole down.
//
// WHAT IS SOLVED AND WHAT IS PRESCRIBED
//
// Honesty about the modelling matters more here than anywhere else in this
// program, because a jet is the one component that a ray tracer cannot get for
// free from the metric:
//
//   * The jet *power* is real physics: it comes out of the Blandford-Znajek
//     formula above, evaluated at the spin and accretion rate of the scene.
//   * The *beaming* is exact.  The observed brightness is computed from the
//     redshift factor g = nu_obs/nu_emit of the actual traced null geodesic
//     contracted with the actual fluid four-velocity, so relativistic
//     boosting, aberration, light bending and gravitational redshift all come
//     out of the same machinery that renders the disc.  The one-sidedness of
//     the jet is therefore a *result*, not an input.
//   * The *radiative transfer* is exact for an optically thin, isotropically
//     emitting plasma: the specific intensity accumulates as
//     int j_nu0 g^(2+alpha) dlambda along the geodesic (derived below).
//   * The jet *geometry, velocity field and emissivity profile* are
//     prescribed.  Solving them would need a full GRMHD simulation, which is
//     a different program.  Each prescription is taken from a measurement or
//     from a conservation law, and each is named where it is used.
//
// THE TRANSFER INTEGRAL
//
// I_nu/nu^3 and j_nu/nu^2 are Lorentz invariants, so along a null geodesic
//
//     d(I_nu / nu^3) / dlambda  =  j_nu / nu^2 .
//
// For a power-law synchrotron source j_nu = j_0 (nu/nu_0)^-alpha, with the
// affine parameter normalised (as make_photon does) so the camera measures
// unit photon energy, the emitted frequency is nu' = nu_obs/g and the integral
// collapses to
//
//     I_nu(obs)  =  (nu_obs/nu_0)^-alpha  .  Integral j_0 g^(2+alpha) dlambda .
//
// Two things follow.  First, the observed *spectrum* is the same power law as
// the emitted one - a shifted power law is still that power law - so the jet
// has a single fixed colour and only its brightness varies across the image.
// Second, all of the relativistic physics is in the single scalar
// g^(2+alpha), which is what the renderer accumulates.
#pragma once

#include "kerr.h"
#include "constants.h"

#include <array>
#include <cmath>
#include <algorithm>

namespace bh {

// Blandford-Znajek jet power, in watts.
//   mdot_c2 : accreted rest-mass energy per unit time, Mdot c^2, in watts
//   phi     : dimensionless magnetic flux threading the horizon (50 = MAD)
//   kappa   : geometry factor, 0.053 for a split-monopole field
inline double blandford_znajek_power(double a, double mdot_c2,
                                     double phi = 50.0, double kappa = 0.053) {
    const double OmH = horizon_angular_velocity(a);
    const double x   = OmH * OmH;
    const double f   = 1.0 + 1.38 * x - 9.2 * x * x;   // higher-order fit
    return (kappa / (4.0 * M_PI)) * phi * phi * x * std::max(0.0, f) * mdot_c2;
}

// Smooth 0 -> 1 ramp, used to give the jet soft ends.  A hard edge would alias
// badly under the four-point quadrature the ray integral uses.
inline double smoothstep01(double t) {
    t = std::clamp(t, 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

struct Jet {
    bool enabled = false;

    // ---- geometry -----------------------------------------------------------
    // Very-long-baseline interferometry resolves M87's jet from 10 to 10^5
    // gravitational radii and finds a *parabolic* boundary, width ~ z^0.58
    // (Asada & Nakamura 2012), not the conical shape of a ballistic outflow.
    // That is the signature of a magnetically collimated flow, so the same
    // exponent is used here.
    double z_base = 3.0;        // emission starts here, in GM/c^2
    double z_top  = 120.0;      // and is truncated here (a real jet does not stop)
    double R_base = 1.6;        // jet radius at z_base
    double k      = 0.58;       // R(z) = R_base (z/z_base)^k

    // ---- kinematics ---------------------------------------------------------
    // Magnetic acceleration converts Poynting flux into bulk motion gradually,
    // and does so in step with the collimation: the asymptotic result of the
    // theory, and what is measured along M87's jet, is Gamma proportional to
    // the jet radius until the flow saturates (Komissarov et al. 2007;
    // Park et al. 2019).  Gamma therefore follows the same z^k law as R.
    double gamma_max = 10.0;
    double z_sat     = 80.0;    // height at which Gamma reaches gamma_max

    // ---- radiation ----------------------------------------------------------
    double alpha = 0.7;         // synchrotron spectral index, S_nu ~ nu^-alpha
    // Jets are observed to be limb-brightened: the emission comes from a sheath
    // around a faster, darker spine (Kim et al. 2018 for M87).  A Gaussian
    // annulus at u = rho/R reproduces that and has no hard edge.
    double u_peak = 0.82, u_wid = 0.17;
    double spine = 0.06, spine_wid = 0.40;
    // Distributed re-acceleration index (see shape()).  zeta = 0 is a jet whose
    // electrons are injected once at the base and thereafter only expand and
    // cool; that jet fades by six orders of magnitude over the length drawn
    // here, and it is ruled out by the plain fact that real jets stay visible
    // out to 10^5 gravitational radii.  Something keeps re-energising the
    // particles along the flow - internal shocks, magnetic reconnection - and
    // zeta is the one number in this model that stands in for it.
    double zeta = 2.0;

    // ---- set by configure() -------------------------------------------------
    double a = 0.0;             // black hole spin, copied in
    double Omega_F = 0.0;       // field-line angular velocity = Omega_H / 2
    double P_jet = 0.0;         // Blandford-Znajek power, W
    double L_rad = 0.0;         // the radiated fraction of it, W
    double eps_rad = 0.02;      // that fraction
    double phi_flux = 50.0;     // dimensionless horizon flux
    double j_nu0 = 0.0;         // emissivity coefficient x r_g,
                                //   W m^-2 sr^-1 Hz^-1 per unit shape()
    double r_bound = 0.0;       // bounding sphere radius, for culling
    double volume_rg3 = 0.0;    // effective emitting volume, in r_g^3

    // Reference frequency for the power law: 550 nm, the middle of the visible
    // band, so the normalisation is never an extrapolation of the fit.
    static constexpr double nu_0 = phys::c / 550e-9;
    // Synchrotron band over which the jet's radiated power is counted.  Blazar
    // and radio-galaxy synchrotron components span roughly radio to soft X-ray;
    // this is what "the jet radiates L_rad" means quantitatively.
    static constexpr double nu_lo = 1.0e9, nu_hi = 1.0e17;

    // ---- geometry helpers ---------------------------------------------------
    double radius(double z_abs) const {
        return R_base * std::pow(std::max(z_abs, 1e-6) / z_base, k);
    }
    double gamma_at(double z_abs) const {
        const double t = std::min(1.0, std::pow(std::max(z_abs, 1e-6) / z_sat, k));
        return 1.0 + (gamma_max - 1.0) * t;
    }
    double speed_at(double z_abs) const {
        const double G = gamma_at(z_abs);
        return std::sqrt(std::max(0.0, 1.0 - 1.0 / (G * G)));
    }

    // Radial emissivity profile across the jet, u = rho / R(z).
    double radial_profile(double u) const {
        const double s = (u - u_peak) / u_wid;
        const double t = u / spine_wid;
        return std::exp(-s * s) + spine * std::exp(-t * t);
    }

    // Dimensionless comoving emissivity at cylindrical (z, rho), both in GM/c^2.
    //
    // The scaling with height is *derived*, not chosen.  Along a steady jet:
    //   - rest mass conservation gives the comoving number density
    //         n ~ 1 / (Gamma R^2),
    //   - flux freezing of the dominant toroidal field gives
    //         B' ~ 1 / (Gamma R),
    //   - and synchrotron emissivity from a power-law electron distribution is
    //         j_nu ~ n B'^(1+alpha) nu^-alpha.
    // Together:   j ~ Gamma^-(2+alpha) R^-(3+alpha).
    //
    // On its own that law is far too steep - see `zeta` above - so the electron
    // normalisation is allowed to grow as (z/z_base)^zeta to represent in-situ
    // re-acceleration.  With the default it gives a ridge-line brightness
    // falling roughly as z^-1, which is what very-long-baseline
    // interferometry actually measures along M87's jet.
    double shape(double z, double rho) const {
        const double az = std::fabs(z);
        if (az <= z_base || az >= z_top) return 0.0;
        const double R = radius(az);
        const double f = radial_profile(rho / R);
        if (f < 1e-9) return 0.0;
        const double g = std::pow(gamma_at(az) / gamma_at(z_base), -(2.0 + alpha)) *
                         std::pow(R / R_base, -(3.0 + alpha)) *
                         std::pow(az / z_base, zeta);
        const double ramp = smoothstep01((az - z_base) / (0.5 * z_base)) *
                            smoothstep01((z_top - az) / (0.15 * z_top));
        return f * g * ramp;
    }

    // Four-velocity of the jet plasma at this point.  Built in the ZAMO frame,
    // which is regular right down to the horizon (a static frame is not - there
    // is no static observer inside the ergosphere, and the jet is launched from
    // deep inside it).  Because the tetrad is orthonormal and the three
    // components form a unit 3-vector times v, the result satisfies
    // u.u = -1 exactly by construction, not to within a tolerance.
    std::array<double, 4> four_velocity(const Geom& g) const {
        const Tetrad T = zamo_tetrad(g);
        const double lapse = zamo_lapse(g);
        const double omega = zamo_omega(g);
        const double sqrt_gpp = g.sn * std::sqrt(g.A / g.Sigma);

        const double z = g.r * g.cs;
        const double v = speed_at(std::fabs(z));

        // The field lines rotate rigidly at Omega_F, which the Blandford-Znajek
        // solution fixes at half the horizon's own angular velocity.  Inside
        // the light cylinder the plasma corotates with them; beyond it, it
        // cannot, and the flow becomes poloidally dominated.  The cap below is
        // the modelling statement that stands in for solving the wind equation.
        double v_phi = (Omega_F - omega) * sqrt_gpp / std::max(lapse, 1e-12);
        v_phi = std::clamp(v_phi, -0.5 * v, 0.5 * v);
        const double v_pol = std::sqrt(std::max(0.0, v * v - v_phi * v_phi));

        // Outflow, away from the equator, in both hemispheres:
        //   +z_hat = cos(th) r_hat - sin(th) th_hat, and -z_hat below.
        const double sgn = (g.cs >= 0.0) ? 1.0 : -1.0;
        const double d_r  =  sgn * g.cs;      // = |cos(th)|
        const double d_th = -sgn * g.sn;

        const double G = 1.0 / std::sqrt(std::max(1e-12, 1.0 - v * v));
        std::array<double, 4> u{};
        for (int mu = 0; mu < 4; ++mu)
            u[mu] = G * (T.e[0][mu] + v_pol * (d_r * T.e[1][mu] + d_th * T.e[2][mu])
                                    + v_phi * T.e[3][mu]);
        return u;
    }

    // Cylindrical coordinates matching bl_to_cartesian: the pseudo-Cartesian
    // (Kerr-Schild style) embedding of Boyer-Lindquist coordinates.
    static void cylindrical(double a_spin, double r, double th, double& z, double& rho) {
        z   = r * std::cos(th);
        rho = std::sqrt(r * r + a_spin * a_spin) * std::sin(th);
    }

    // -------------------------------------------------------------------------
    // Tie the emissivity to the Blandford-Znajek power.
    //
    //   mdot_c2 : Mdot c^2 of the accretion flow, in watts
    //   r_g_m   : the gravitational radius in metres
    //
    // The emissivity amplitude is fixed by requiring that the comoving
    // emission integrated over both jets equals L_rad = eps_rad * P_jet:
    //
    //     4 pi Integral j_bol dV  =  L_rad .
    // -------------------------------------------------------------------------
    void configure(double spin, double mdot_c2, double r_g_m) {
        a = spin;
        Omega_F = 0.5 * horizon_angular_velocity(spin);
        P_jet = blandford_znajek_power(spin, mdot_c2, phi_flux);
        L_rad = eps_rad * P_jet;
        r_bound = std::sqrt(z_top * z_top + radius(z_top) * radius(z_top)) * 1.05;

        // Transverse integral of the radial profile, Integral f(u) u du.
        double C_f = 0.0;
        {
            constexpr int N = 2000;
            const double du = 3.0 / N;
            for (int i = 0; i < N; ++i) {
                const double u = (i + 0.5) * du;
                C_f += radial_profile(u) * u * du;
            }
        }
        // Longitudinal integral, with the transverse one folded in.  Both jets.
        double V = 0.0;
        {
            constexpr int N = 4000;
            const double dz = (z_top - z_base) / N;
            for (int i = 0; i < N; ++i) {
                const double z = z_base + (i + 0.5) * dz;
                const double R = radius(z);
                // shape at rho = 0 would sit in the dark spine, so take the
                // z-dependent part directly.
                const double gz = std::pow(gamma_at(z) / gamma_at(z_base), -(2.0 + alpha)) *
                                  std::pow(R / R_base, -(3.0 + alpha)) *
                                  std::pow(z / z_base, zeta);
                const double ramp = smoothstep01((z - z_base) / (0.5 * z_base)) *
                                    smoothstep01((z_top - z) / (0.15 * z_top));
                V += gz * ramp * R * R * dz;
            }
            V *= 2.0 * (2.0 * M_PI) * C_f;   // 2 pi from phi, factor 2 for both jets
        }
        volume_rg3 = V;
        const double V_m3 = V * r_g_m * r_g_m * r_g_m;
        if (!(V_m3 > 0.0) || !(L_rad > 0.0)) { j_nu0 = 0.0; return; }

        // Bolometric emissivity amplitude, W m^-3 sr^-1.
        const double A = L_rad / (4.0 * M_PI * V_m3);

        // Convert it to the coefficient of the power law at nu_0:
        //   j_bol = j_nu0 nu_0^alpha Integral_{nu_lo}^{nu_hi} nu^-alpha dnu.
        const double p = 1.0 - alpha;
        const double band = (std::fabs(p) > 1e-9)
            ? (std::pow(nu_hi, p) - std::pow(nu_lo, p)) / p
            : std::log(nu_hi / nu_lo);
        j_nu0 = A / (std::pow(nu_0, alpha) * band) * r_g_m;
    }

    // Observed intensity is proportional to g^(2+alpha); see the header note.
    double boost(double g) const { return std::pow(g, 2.0 + alpha); }
};

} // namespace bh
