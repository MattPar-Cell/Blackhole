// spectrum.h - Physical blackbody radiance -> CIE XYZ -> linear sRGB.
//
// Nothing here is an artistic colour ramp.  A temperature is turned into a
// colour by integrating the Planck spectral radiance against the CIE 1931
// 2-degree standard observer, so the disc, the Sun and the stars all end up on
// one absolute radiometric scale (W m^-2 sr^-1) and their relative brightness
// in the final image is the real one.
//
// The key relativistic identity used by the renderer is that the photon
// occupation number I_nu / nu^3 is a Lorentz invariant conserved along a null
// geodesic (Liouville's theorem).  A blackbody of temperature T_emit therefore
// arrives at the camera as an *exact blackbody* of temperature
//
//     T_obs = g * T_emit,      g = (p.u)_obs / (p.u)_emit
//
// with no extra factor: the familiar g^4 boost of the bolometric intensity is
// already contained in the sigma T^4 scaling of the Planck function.  So the
// renderer only ever has to ask for the colour of a blackbody at g*T.
#pragma once

namespace spec {

struct XYZ {
    double x = 0.0, y = 0.0, z = 0.0;
};

struct RGB {
    double r = 0.0, g = 0.0, b = 0.0;
};

// Planck spectral radiance B_lambda(T), in W m^-2 sr^-1 m^-1.
double planck_radiance(double lambda_m, double T);

// Integrate B_lambda(T) against the CIE 1931 colour matching functions over
// 360-830 nm.  Result is in W m^-2 sr^-1 (radiometrically scaled, not
// normalised), so magnitudes are physically comparable between sources.
XYZ blackbody_xyz_exact(double T);

// Table-accelerated version of the above (log-log interpolation, <0.1% error).
// Safe for T from 50 K to 3e9 K; clamps outside that range.
XYZ blackbody_xyz(double T);

// Total bolometric radiance of a blackbody, sigma T^4 / pi  (W m^-2 sr^-1).
double blackbody_bolometric(double T);

// CIE XYZ (D65) -> linear sRGB primaries.  May produce out-of-gamut negatives.
RGB xyz_to_linear_srgb(const XYZ& c);

// Desaturate towards D65 white until every channel is non-negative, keeping
// luminance fixed.  This is what makes deep-red and deep-violet blackbodies
// display sensibly instead of clipping to black.
RGB gamut_clamp(const RGB& c);

// Convenience: physical radiance of a blackbody at T, as linear sRGB.
RGB blackbody_rgb(double T);

// Approximate effective temperature of a main-sequence star from its B-V
// colour index (Ballesteros 2012).
double bv_to_temperature(double bv);

} // namespace spec
