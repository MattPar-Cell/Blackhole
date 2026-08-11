// eos.h - Equations of state for degenerate matter.
//
// The structure of a neutron star is set by one relation: how pressure
// responds to density, P(rho).  That is the whole physical input to the TOV
// equation, and it is where all the uncertainty in neutron star physics lives.
//
// Three of the four equations of state here are derived from first principles,
// which means they can be validated against exact results rather than trusted:
//
//   IdealFermiGas   - a free degenerate gas, from Fermi-Dirac statistics at
//                     T = 0.  With neutrons this reproduces the historical
//                     Oppenheimer-Volkoff limit of 0.71 solar masses; with
//                     electrons it reproduces Chandrasekhar's 1.44.
//   Polytrope       - P = K rho^Gamma, which in the Newtonian limit has an
//                     exact Lane-Emden solution to check against.
//   UniformDensity  - incompressible matter, whose relativistic structure is
//                     the exact interior Schwarzschild solution.
//
// The fourth is a piecewise polytrope, the standard way of representing a
// realistic tabulated equation of state with a handful of numbers (Read,
// Lackey, Owen & Friedman 2009).  Its parameters are fits to published
// microphysical calculations rather than anything derived here, and the
// honest test of those is whether the resulting stars match what pulsar
// timing and NICER actually measure.
#pragma once

#include "constants.h"

#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

namespace ns {

// An equation of state supplies pressure from mass-energy density and back.
// Both are in SI: rho in kg m^-3 (total mass-energy density divided by c^2),
// P in pascals.
class EOS {
public:
    virtual ~EOS() = default;
    virtual double pressure(double rho) const = 0;   // P(rho)
    virtual double density(double P) const = 0;      // rho(P), the inverse
    virtual std::string name() const = 0;

    // Speed of sound squared, in units of c^2:  c_s^2 = dP/d(rho c^2).
    // Causality requires this to stay below 1; the TOV solver checks it.
    virtual double sound_speed_squared(double rho) const {
        const double h = 1e-4 * rho;
        return (pressure(rho + h) - pressure(rho - h)) / (2.0 * h * phys::c * phys::c);
    }
};

// ---------------------------------------------------------------------------
// Ideal degenerate Fermi gas at zero temperature.
//
// For a gas of fermions of mass m filling momentum states up to p_F, write
// x = p_F / (m c).  Integrating over the filled Fermi sphere gives
//
//   n   = (8 pi / 3) (m c / h)^3 x^3
//   eps = K [ (2x^3 + x) sqrt(1+x^2) - asinh(x) ]
//   P   = (K/3) [ (2x^3 - 3x) sqrt(1+x^2) + 3 asinh(x) ]
//
// with K = pi m^4 c^5 / h^3.  These satisfy dE/dn = m c^2 sqrt(1+x^2) and
// P = mu n - eps exactly, which validate() checks.
//
// `mass_per_fermion` separates the two roles mass plays: the pressure comes
// from the degenerate species, but the inertia can come from something else.
// In a white dwarf the electrons hold the star up while the nucleons carry
// essentially all the mass, so mass_per_fermion = mu_e * m_u with mu_e = 2.
// ---------------------------------------------------------------------------
class IdealFermiGas : public EOS {
public:
    IdealFermiGas(double m_fermion, double mass_per_fermion, std::string label)
        : m_(m_fermion), mu_(mass_per_fermion), label_(std::move(label)) {
        K_ = M_PI * m_ * m_ * m_ * m_ * std::pow(phys::c, 5) /
             (phys::h_planck * phys::h_planck * phys::h_planck);
        n_scale_ = (8.0 * M_PI / 3.0) *
                   std::pow(m_ * phys::c / phys::h_planck, 3.0);
        // Whether the fermions' own kinetic energy counts towards the inertia:
        // it does for neutrons (they are the mass), not meaningfully for the
        // electrons in a white dwarf.
        self_gravitating_ = std::fabs(mu_ - m_) < 1e-3 * m_;
    }

    double number_density(double x) const { return n_scale_ * x * x * x; }

    double pressure_of_x(double x) const {
        const double s = std::sqrt(1.0 + x * x);
        return (K_ / 3.0) * ((2.0 * x * x * x - 3.0 * x) * s + 3.0 * std::asinh(x));
    }

    // Total mass-energy density divided by c^2.
    double density_of_x(double x) const {
        if (self_gravitating_) {
            const double s = std::sqrt(1.0 + x * x);
            const double eps = K_ * ((2.0 * x * x * x + x) * s - std::asinh(x));
            return eps / (phys::c * phys::c);
        }
        return mu_ * number_density(x);
    }

    double chemical_potential(double x) const {
        return m_ * phys::c * phys::c * std::sqrt(1.0 + x * x);
    }

    double pressure(double rho) const override {
        return pressure_of_x(x_of_density(rho));
    }

    double density(double P) const override {
        if (P <= 0.0) return 0.0;
        // Bracket then bisect: P(x) is strictly increasing.
        double lo = 1e-12, hi = 1.0;
        while (pressure_of_x(hi) < P && hi < 1e12) hi *= 2.0;
        for (int i = 0; i < 200; ++i) {
            const double mid = std::sqrt(lo * hi);      // geometric bisection
            if (pressure_of_x(mid) < P) lo = mid; else hi = mid;
        }
        return density_of_x(std::sqrt(lo * hi));
    }

    double x_of_density(double rho) const {
        if (rho <= 0.0) return 0.0;
        double lo = 1e-12, hi = 1.0;
        while (density_of_x(hi) < rho && hi < 1e12) hi *= 2.0;
        for (int i = 0; i < 200; ++i) {
            const double mid = std::sqrt(lo * hi);
            if (density_of_x(mid) < rho) lo = mid; else hi = mid;
        }
        return std::sqrt(lo * hi);
    }

    std::string name() const override { return label_; }

private:
    double m_, mu_, K_, n_scale_;
    bool self_gravitating_;
    std::string label_;
};

inline IdealFermiGas neutron_gas() {
    return IdealFermiGas(phys::m_n, phys::m_n, "ideal degenerate neutron gas");
}
// mu_e = 2 for carbon/oxygen: two nucleons per electron.
inline IdealFermiGas electron_gas(double mu_e = 2.0) {
    return IdealFermiGas(phys::m_e, mu_e * phys::m_u, "ideal degenerate electron gas");
}

// ---------------------------------------------------------------------------
// P = K rho^Gamma.
// ---------------------------------------------------------------------------
class Polytrope : public EOS {
public:
    Polytrope(double K, double Gamma) : K_(K), G_(Gamma) {}
    double pressure(double rho) const override {
        return (rho > 0.0) ? K_ * std::pow(rho, G_) : 0.0;
    }
    double density(double P) const override {
        return (P > 0.0) ? std::pow(P / K_, 1.0 / G_) : 0.0;
    }
    std::string name() const override {
        return "polytrope Gamma = " + std::to_string(G_);
    }
private:
    double K_, G_;
};

// ---------------------------------------------------------------------------
// Incompressible matter: rho = const wherever P > 0.  Physically absurd (the
// sound speed is infinite) but it is the one case where the TOV equation has a
// closed-form solution, the interior Schwarzschild metric of 1916, so it is
// the sharpest possible test of the integrator.
// ---------------------------------------------------------------------------
class UniformDensity : public EOS {
public:
    // Incompressible matter has no P(rho) to invert - the density is the same
    // everywhere the pressure is positive - so the central pressure has to be
    // supplied.  Passing the analytic value for a chosen radius turns the TOV
    // integration into a direct test: it should come back out at that radius.
    UniformDensity(double rho, double P_central) : rho_(rho), Pc_(P_central) {}
    double pressure(double rho) const override { return (rho >= rho_) ? Pc_ : 0.0; }
    double density(double P) const override { return (P > 0.0) ? rho_ : 0.0; }
    double sound_speed_squared(double) const override { return 1e30; }
    std::string name() const override { return "uniform density"; }
    double rho() const { return rho_; }
private:
    double rho_, Pc_;
};

// ---------------------------------------------------------------------------
// Piecewise polytrope (Read, Lackey, Owen & Friedman 2009).
//
// A realistic tabulated equation of state is reproduced to a few percent by
// three polytropic segments above a fixed crust, parameterised by the pressure
// p1 at rho = 10^14.7 g/cm^3 and the three adiabatic indices.  The segments are
// joined so that P is continuous; the density-energy relation is taken in the
// same polytropic form, which is the standard simplification.
//
// The parameter sets below are literature fits, not derived here.  What can be
// checked - and validate() does check it - is that the resulting stars have the
// masses and radii that pulsar timing and NICER measure.
// ---------------------------------------------------------------------------
class PiecewisePolytrope : public EOS {
public:
    PiecewisePolytrope(double log_p1_cgs, double G1, double G2, double G3,
                       std::string label)
        : label_(std::move(label)) {
        const double c2 = phys::c * phys::c;

        // Dividing *rest-mass* densities, converted from g/cm^3 to kg/m^3.
        rho_[1] = std::pow(10.0, 14.7) * 1000.0;
        rho_[2] = std::pow(10.0, 15.0) * 1000.0;
        gamma_[1] = G1; gamma_[2] = G2; gamma_[3] = G3;

        // p1 is quoted in dyn/cm^2; 1 dyn/cm^2 = 0.1 Pa.
        const double p1 = std::pow(10.0, log_p1_cgs) * 0.1;
        K_[1] = p1 / std::pow(rho_[1], gamma_[1]);

        // Continuity of P at each dividing density fixes the other constants.
        K_[2] = K_[1] * std::pow(rho_[1], gamma_[1] - gamma_[2]);
        K_[3] = K_[2] * std::pow(rho_[2], gamma_[2] - gamma_[3]);

        // Crust: a single polytrope, P = K rho^1.35692, with K = 3.99874e-8 in
        // cgs (Read et al. 2009).  Converting, P_SI = 0.1 K_cgs (rho_SI/1000)^G.
        gamma_[0] = 1.35692;
        K_[0] = 0.1 * 3.99873692e-8 * std::pow(1000.0, -gamma_[0]);
        // The crust meets the first core segment where the two pressures agree.
        rho_[0] = std::pow(K_[0] / K_[1], 1.0 / (gamma_[1] - gamma_[0]));

        // Energy density in this parametrisation is
        //     eps = (1 + a_i) rho c^2 + P / (Gamma_i - 1)
        // where rho is *rest-mass* density and P/(Gamma-1) is the internal
        // energy.  Skipping that internal-energy term is not a small error:
        // it under-counts what gravitates, which makes the star behave as if
        // its equation of state were stiffer and inflates the maximum mass by
        // about 20%.  The constants a_i follow from requiring eps itself to be
        // continuous where the segments meet.
        a_[0] = 0.0;
        for (int i = 0; i < 3; ++i) {
            a_[i + 1] = a_[i]
                + K_[i]     * std::pow(rho_[i], gamma_[i] - 1.0)     / ((gamma_[i] - 1.0) * c2)
                - K_[i + 1] * std::pow(rho_[i], gamma_[i + 1] - 1.0) / ((gamma_[i + 1] - 1.0) * c2);
        }
    }

    // P from rest-mass density.
    double pressure_of_rest(double rho_rest) const {
        if (rho_rest <= 0.0) return 0.0;
        const int i = segment_of_rest(rho_rest);
        return K_[i] * std::pow(rho_rest, gamma_[i]);
    }

    // Total mass-energy density / c^2, from rest-mass density.
    double energy_of_rest(double rho_rest) const {
        if (rho_rest <= 0.0) return 0.0;
        const int i = segment_of_rest(rho_rest);
        const double P = K_[i] * std::pow(rho_rest, gamma_[i]);
        return (1.0 + a_[i]) * rho_rest + P / ((gamma_[i] - 1.0) * phys::c * phys::c);
    }

    double density(double P) const override {
        if (P <= 0.0) return 0.0;
        return energy_of_rest(rest_of_pressure(P));
    }

    double pressure(double rho_energy) const override {
        if (rho_energy <= 0.0) return 0.0;
        // Invert eps(rho_rest), which is monotonic.
        double lo = 1e-6 * rho_energy, hi = rho_energy;
        while (energy_of_rest(hi) < rho_energy && hi < 1e25) hi *= 2.0;
        for (int k = 0; k < 200; ++k) {
            const double mid = std::sqrt(lo * hi);
            if (energy_of_rest(mid) < rho_energy) lo = mid; else hi = mid;
        }
        return pressure_of_rest(std::sqrt(lo * hi));
    }

    std::string name() const override { return label_; }

private:
    int segment_of_rest(double rho) const {
        if (rho < rho_[0]) return 0;
        if (rho < rho_[1]) return 1;
        if (rho < rho_[2]) return 2;
        return 3;
    }
    double rest_of_pressure(double P) const {
        for (int i = 0; i < 4; ++i) {
            const double rho = std::pow(P / K_[i], 1.0 / gamma_[i]);
            if (segment_of_rest(rho) == i) return rho;
        }
        return std::pow(P / K_[3], 1.0 / gamma_[3]);
    }
    double K_[4], gamma_[4], a_[4], rho_[3];
    std::string label_;
};

// Two literature parameter sets spanning the plausible range: SLy is a
// moderately soft equation of state that has been a standard reference for
// decades, MS1 a stiff one.  A stiffer equation of state supports more mass
// and gives larger stars.
inline PiecewisePolytrope eos_sly() {
    return PiecewisePolytrope(34.384, 3.005, 2.988, 2.851, "SLy (piecewise polytrope fit)");
}
inline PiecewisePolytrope eos_stiff() {
    return PiecewisePolytrope(34.858, 3.224, 3.033, 1.325, "MS1 (stiff, piecewise polytrope fit)");
}

} // namespace ns
