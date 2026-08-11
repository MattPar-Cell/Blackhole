// constants.h - CODATA 2018 / IAU 2015 physical constants (SI) and unit helpers.
//
// The simulation integrates geodesics in *geometrised units* where
//     G = c = M = 1,
// so that lengths, times and masses are all measured in units of the
// gravitational radius  r_g = GM/c^2.  Everything that has to be reported to a
// human (temperatures, luminosities, angular sizes, the size of the Sun) is
// converted back to SI with the helpers at the bottom of this file.
#pragma once

#include <cmath>

// M_PI is a POSIX extension, not standard C++.  MSVC only provides it when
// _USE_MATH_DEFINES is set before <cmath>, so define it here if it is missing.
// Every translation unit that uses M_PI includes this header.
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace phys {

// ---- Fundamental constants (SI) --------------------------------------------
inline constexpr double G          = 6.67430e-11;      // m^3 kg^-1 s^-2
inline constexpr double c          = 2.99792458e8;     // m s^-1  (exact)
inline constexpr double h_planck   = 6.62607015e-34;   // J s     (exact)
inline constexpr double k_B        = 1.380649e-23;     // J K^-1  (exact)
inline constexpr double sigma_SB   = 5.670374419e-8;   // W m^-2 K^-4
inline constexpr double sigma_T    = 6.6524587321e-29; // m^2, Thomson cross-section
inline constexpr double m_p        = 1.67262192369e-27;// kg, proton
inline constexpr double m_n        = 1.67492749804e-27;// kg, neutron
inline constexpr double m_e        = 9.1093837015e-31; // kg, electron
inline constexpr double m_u        = 1.66053906660e-27;// kg, atomic mass unit
inline constexpr double eV         = 1.602176634e-19;  // J (exact)
inline constexpr double MeV        = 1.0e6 * eV;

// Nuclear saturation density: the density of matter inside a heavy nucleus,
// about 0.16 baryons per cubic femtometre.  Neutron star cores run at a few
// times this.
inline constexpr double rho_nuclear = 2.7e17;          // kg m^-3

// ---- Solar / astronomical (IAU 2015 nominal values) ------------------------
inline constexpr double M_sun      = 1.98847e30;       // kg
inline constexpr double R_sun      = 6.957e8;          // m   (nominal solar radius)
inline constexpr double L_sun      = 3.828e26;         // W   (nominal solar luminosity)
inline constexpr double T_sun      = 5772.0;           // K   (solar effective temperature)
inline constexpr double AU         = 1.495978707e11;   // m   (exact)
inline constexpr double parsec     = 3.0856775814913673e16; // m
inline constexpr double lightyear  = 9.4607304725808e15;    // m
inline constexpr double year       = 3.15576e7;        // s (Julian year)

// Gravitational radius r_g = GM/c^2 and Schwarzschild radius r_s = 2 r_g.
inline double r_g_metres(double M_kg)  { return G * M_kg / (c * c); }
inline double r_s_metres(double M_kg)  { return 2.0 * G * M_kg / (c * c); }

// Eddington luminosity for fully ionised hydrogen: L_Edd = 4 pi G M m_p c / sigma_T.
inline double L_eddington(double M_kg) {
    return 4.0 * M_PI * G * M_kg * m_p * c / sigma_T;
}

// Eddington accretion rate for a disc of radiative efficiency eta:
//   Mdot_Edd = L_Edd / (eta c^2)
inline double mdot_eddington(double M_kg, double eta) {
    return L_eddington(M_kg) / (eta * c * c);
}

// Angular diameter (radians) of an object of radius R seen from distance D,
// ignoring lensing:  theta = 2 arcsin(R/D).
inline double angular_diameter(double R, double D) {
    if (D <= R) return M_PI;
    return 2.0 * std::asin(R / D);
}

inline double rad_to_arcsec(double x)      { return x * 180.0 / M_PI * 3600.0; }
inline double rad_to_microarcsec(double x) { return rad_to_arcsec(x) * 1.0e6; }

} // namespace phys
