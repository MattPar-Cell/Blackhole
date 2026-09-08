#include "spectrum.h"
#include "constants.h"

#include <cmath>
#include <algorithm>
#include <vector>
#include <mutex>

namespace spec {

using namespace phys;

double planck_radiance(double lambda_m, double T) {
    if (T <= 0.0 || lambda_m <= 0.0) return 0.0;
    const double l  = lambda_m;
    const double l5 = l * l * l * l * l;
    const double x  = h_planck * c / (l * k_B * T);
    // expm1 keeps the Rayleigh-Jeans limit (x << 1) accurate; for x > 700 the
    // Wien tail is zero to double precision.
    if (x > 700.0) return 0.0;
    return (2.0 * h_planck * c * c) / (l5 * std::expm1(x));
}

double blackbody_bolometric(double T) {
    return sigma_SB * T * T * T * T / M_PI;
}

// ---------------------------------------------------------------------------
// CIE 1931 2-degree colour matching functions.
//
// Multi-lobe piecewise-Gaussian fit of Wyman, Sloan & Shirley (2013),
// "Simple Analytic Approximations to the CIE XYZ Color Matching Functions",
// Journal of Computer Graphics Techniques 2(2).  Maximum error under 1% of
// peak, which is far below the accuracy of anything else in this renderer.
// ---------------------------------------------------------------------------
static inline double lobe(double x, double mu, double s1, double s2) {
    const double t = (x - mu) * ((x < mu) ? 1.0 / s1 : 1.0 / s2);
    return std::exp(-0.5 * t * t);
}

static inline double cie_x(double nm) {
    return 1.056 * lobe(nm, 599.8, 37.9, 31.0)
         + 0.362 * lobe(nm, 442.0, 16.0, 26.7)
         - 0.065 * lobe(nm, 501.1, 20.4, 26.2);
}
static inline double cie_y(double nm) {
    return 0.821 * lobe(nm, 568.8, 46.9, 40.5)
         + 0.286 * lobe(nm, 530.9, 16.3, 31.1);
}
static inline double cie_z(double nm) {
    return 1.217 * lobe(nm, 437.0, 11.8, 36.0)
         + 0.681 * lobe(nm, 459.0, 26.0, 13.8);
}

XYZ blackbody_xyz_exact(double T) {
    XYZ out;
    if (T <= 0.0) return out;
    // 1 nm steps across the full range where the CMFs are non-negligible.
    // The integration variable is wavelength in metres so that the result
    // inherits the units of B_lambda, i.e. W m^-2 sr^-1: Y is the visible-band
    // radiance weighted by the eye's photopic response, and X:Y:Z is the
    // chromaticity.  No further normalisation is applied - the one overall
    // constant is common to every source in the scene and is absorbed by the
    // exposure control, so relative brightnesses stay physical.
    constexpr double lo = 360.0, hi = 830.0, dl = 1.0;
    for (double nm = lo; nm <= hi; nm += dl) {
        const double w = planck_radiance(nm * 1e-9, T) * dl * 1e-9;
        out.x += w * cie_x(nm);
        out.y += w * cie_y(nm);
        out.z += w * cie_z(nm);
    }
    return out;
}

XYZ cie_cmf(double nm) { return {cie_x(nm), cie_y(nm), cie_z(nm)}; }

XYZ powerlaw_xyz(double alpha) {
    // I_nu = A (nu/nu_0)^-alpha  ->  I_lambda = I_nu c / lambda^2
    //                                        = A c (lambda/lambda_0)^alpha / lambda^2
    constexpr double lambda_0 = 550e-9;
    XYZ out;
    constexpr double lo = 360.0, hi = 830.0, dl = 1.0;
    for (double nm = lo; nm <= hi; nm += dl) {
        const double lam = nm * 1e-9;
        const double I_l = c * std::pow(lam / lambda_0, alpha) / (lam * lam);
        const double w = I_l * dl * 1e-9;
        out.x += w * cie_x(nm);
        out.y += w * cie_y(nm);
        out.z += w * cie_z(nm);
    }
    return out;
}

// --- log-log lookup table --------------------------------------------------
namespace {

constexpr int    kTableN    = 4096;
constexpr double kLogTMin   = 1.6989700043360187;   // log10(50 K)
constexpr double kLogTMax   = 9.4771212547196626;   // log10(3e9 K)

struct Table {
    std::vector<double> lx, ly, lz;   // natural logs
    Table() {
        lx.resize(kTableN); ly.resize(kTableN); lz.resize(kTableN);
        for (int i = 0; i < kTableN; ++i) {
            const double f = static_cast<double>(i) / (kTableN - 1);
            const double T = std::pow(10.0, kLogTMin + f * (kLogTMax - kLogTMin));
            const XYZ v = blackbody_xyz_exact(T);
            lx[i] = std::log(std::max(v.x, 1e-300));
            ly[i] = std::log(std::max(v.y, 1e-300));
            lz[i] = std::log(std::max(v.z, 1e-300));
        }
    }
};

const Table& table() {
    static const Table t;   // thread-safe static init (C++11 magic statics)
    return t;
}

} // namespace

XYZ blackbody_xyz(double T) {
    const Table& t = table();
    const double lg = std::log10(std::clamp(T, 50.0, 3.0e9));
    const double f  = (lg - kLogTMin) / (kLogTMax - kLogTMin) * (kTableN - 1);
    const int    i  = std::clamp(static_cast<int>(f), 0, kTableN - 2);
    const double u  = f - i;
    XYZ o;
    o.x = std::exp(t.lx[i] * (1.0 - u) + t.lx[i + 1] * u);
    o.y = std::exp(t.ly[i] * (1.0 - u) + t.ly[i + 1] * u);
    o.z = std::exp(t.lz[i] * (1.0 - u) + t.lz[i + 1] * u);
    return o;
}

// sRGB / Rec.709 primaries with a D65 white point (IEC 61966-2-1).
RGB xyz_to_linear_srgb(const XYZ& c) {
    RGB o;
    o.r =  3.2404542 * c.x - 1.5371385 * c.y - 0.4985314 * c.z;
    o.g = -0.9692660 * c.x + 1.8760108 * c.y + 0.0415560 * c.z;
    o.b =  0.0556434 * c.x - 0.2040259 * c.y + 1.0572252 * c.z;
    return o;
}

RGB gamut_clamp(const RGB& c) {
    const double m = std::min({c.r, c.g, c.b});
    if (m >= 0.0) return c;
    // Adding equal amounts of all three primaries moves the colour along the
    // line towards white without changing its position on the luminance axis
    // by more than the (small) amount that was out of gamut to begin with.
    return {c.r - m, c.g - m, c.b - m};
}

RGB blackbody_rgb(double T) {
    return gamut_clamp(xyz_to_linear_srgb(blackbody_xyz(T)));
}

double bv_to_temperature(double bv) {
    // Ballesteros (2012), EPL 97, 34008.
    const double x = std::clamp(bv, -0.4, 2.0);
    return 4600.0 * (1.0 / (0.92 * x + 1.7) + 1.0 / (0.92 * x + 0.62));
}

} // namespace spec
