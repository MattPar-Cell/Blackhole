// scene.h - The things light can come from: stars, the Sun, and the camera
//           that collects them.
//
// Every source reports a spectral radiance in the same absolute units used by
// spectrum.h (CIE-weighted W m^-2 sr^-1), so the brightness ratios between the
// accretion disc, the solar photosphere and a sixth-magnitude star in the
// final image are the physically correct ones.  Only the display transform in
// image.h compresses that range for a screen.
#pragma once

#include "kerr.h"
#include "geodesic.h"
#include "spectrum.h"
#include "constants.h"

#include <vector>
#include <cmath>
#include <cstdint>
#include <algorithm>

namespace bh {

// ---------------------------------------------------------------------------
// Small deterministic PRNG (splitmix64 + xoshiro-style output) so that scenes
// are reproducible across platforms and thread counts.
// ---------------------------------------------------------------------------
struct Rng {
    uint64_t s;
    explicit Rng(uint64_t seed = 0x9E3779B97F4A7C15ull) : s(seed) {}
    uint64_t next_u64() {
        s += 0x9E3779B97F4A7C15ull;
        uint64_t z = s;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        return z ^ (z >> 31);
    }
    double uniform() { return (next_u64() >> 11) * 0x1.0p-53; }
    double uniform(double a, double b) { return a + (b - a) * uniform(); }
    double normal() {
        // Box-Muller; only one of the pair is kept (cheap enough here).
        const double u1 = std::max(uniform(), 1e-15), u2 = uniform();
        return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
    }
};

// ---------------------------------------------------------------------------
// Star field.
//
// Stars are point sources: each carries a V-band irradiance at the camera
// (W m^-2) derived from an apparent magnitude, and a colour temperature from a
// B-V index.  The magnitude distribution follows the observed cumulative star
// count log N(<m) ~ 0.45 m, which reproduces the real naked-eye sky density
// (about 9000 stars brighter than V = 6.5 over the whole sphere).
//
// Because the stars live at infinity and are looked up with the *final*
// direction of each traced geodesic, they are lensed exactly: the field shows
// the Einstein ring, the secondary images just outside the shadow, and the
// characteristic swirl of the frame-dragged photon ring for free.
// ---------------------------------------------------------------------------
class StarField {
public:
    StarField(int count, uint64_t seed, bool milky_way = true)
        : n_lat_(360), n_lon_(720), milky_way_(milky_way) {
        Rng rng(seed);
        bins_.resize(static_cast<size_t>(n_lat_) * n_lon_);

        stars_.reserve(count);
        for (int i = 0; i < count; ++i) {
            Star s;

            // Direction: isotropic, then a fraction pulled towards the
            // galactic plane to make the Milky Way a real overdensity.
            double dir[3];
            if (rng.uniform() < 0.45) {
                const double b = rng.normal() * 0.12;            // galactic latitude
                const double l = rng.uniform(0.0, 2.0 * M_PI);
                const double cb = std::cos(b);
                double gd[3] = {cb * std::cos(l), cb * std::sin(l), std::sin(b)};
                galactic_to_world(gd, dir);
            } else {
                const double z  = rng.uniform(-1.0, 1.0);
                const double ph = rng.uniform(0.0, 2.0 * M_PI);
                const double rr = std::sqrt(std::max(0.0, 1.0 - z * z));
                dir[0] = rr * std::cos(ph); dir[1] = rr * std::sin(ph); dir[2] = z;
            }
            s.x = dir[0]; s.y = dir[1]; s.z = dir[2];

            // Apparent magnitude by inverse-CDF sampling of N(<m) ~ 10^{k m}.
            constexpr double k = 0.45, m_lo = -1.5, m_hi = 12.0;
            const double u = rng.uniform();
            s.mag = m_lo + std::log10(u * (std::pow(10.0, k * (m_hi - m_lo)) - 1.0) + 1.0) / k;

            // B-V index: a rough but realistic mix of spectral types for
            // magnitude-limited samples (blue-white dominated, red tail).
            double bv = 0.55 + 0.42 * rng.normal();
            if (rng.uniform() < 0.10) bv = rng.uniform(1.2, 1.9);   // red giants
            s.temperature = spec::bv_to_temperature(std::clamp(bv, -0.35, 1.95));

            // V = 0 corresponds to 3.64e-9 W m^-2 um^-1 over an 0.089 um band.
            constexpr double F0 = 3.64e-9 * 0.089 * 1e6 * 1e-6;   // = 3.24e-10 W m^-2
            s.irradiance = F0 * std::pow(10.0, -0.4 * s.mag);

            stars_.push_back(s);
        }

        for (size_t i = 0; i < stars_.size(); ++i) bins_[bin_of(stars_[i])].push_back(static_cast<uint32_t>(i));
    }

    // Radiance seen looking along unit direction d, with a point-spread
    // function of angular width `sigma` (radians, ~0.5 pixel).
    spec::XYZ radiance(const double d[3], double sigma) const {
        spec::XYZ out;

        const double cutoff = 3.0 * sigma;
        const double inv2s2 = 1.0 / (2.0 * sigma * sigma);
        // A normalised 2D Gaussian on the sky: integral over solid angle = 1,
        // so irradiance (W m^-2) becomes radiance (W m^-2 sr^-1).
        const double norm = 1.0 / (2.0 * M_PI * sigma * sigma);

        const double theta = std::acos(std::clamp(d[2], -1.0, 1.0));
        const double phi   = std::atan2(d[1], d[0]);
        const int    it    = static_cast<int>(theta / M_PI * n_lat_);
        const double dlat  = M_PI / n_lat_;
        const int    span_t = std::max(1, static_cast<int>(std::ceil(cutoff / dlat)));
        const double sn    = std::max(std::sin(theta), 1e-6);
        const double dlon  = 2.0 * M_PI / n_lon_;
        const int    span_p = std::max(1, static_cast<int>(std::ceil(cutoff / (sn * dlon))));

        for (int a = -span_t; a <= span_t; ++a) {
            const int ti = it + a;
            if (ti < 0 || ti >= n_lat_) continue;
            // Near the poles the longitude bins are narrow; widen or wrap fully.
            const int sp = (span_p * 2 + 1 >= n_lon_) ? n_lon_ / 2 : span_p;
            for (int b = -sp; b <= sp; ++b) {
                int pi = static_cast<int>((phi + M_PI) / (2.0 * M_PI) * n_lon_) + b;
                pi = ((pi % n_lon_) + n_lon_) % n_lon_;
                for (uint32_t idx : bins_[static_cast<size_t>(ti) * n_lon_ + pi]) {
                    const Star& s = stars_[idx];
                    const double dot = std::clamp(d[0] * s.x + d[1] * s.y + d[2] * s.z, -1.0, 1.0);
                    const double ang = (dot > 0.9999995) ? std::sqrt(std::max(0.0, 2.0 * (1.0 - dot)))
                                                         : std::acos(dot);
                    if (ang > cutoff) continue;
                    const double w = norm * std::exp(-ang * ang * inv2s2) * s.irradiance;
                    const spec::XYZ c = spec::blackbody_xyz(s.temperature);
                    const double inv = (c.y > 0.0) ? 1.0 / c.y : 0.0;
                    out.x += w * c.x * inv;
                    out.y += w;
                    out.z += w * c.z * inv;
                }
            }
        }

        if (milky_way_) {
            const spec::XYZ mw = milky_way_radiance(d);
            out.x += mw.x; out.y += mw.y; out.z += mw.z;
        }
        return out;
    }

    size_t size() const { return stars_.size(); }

private:
    struct Star {
        double x, y, z;
        double mag, temperature, irradiance;
    };

    std::vector<Star> stars_;
    std::vector<std::vector<uint32_t>> bins_;
    int  n_lat_, n_lon_;
    bool milky_way_;

    size_t bin_of(const Star& s) const {
        const double theta = std::acos(std::clamp(s.z, -1.0, 1.0));
        const double phi   = std::atan2(s.y, s.x);
        int ti = static_cast<int>(theta / M_PI * n_lat_);
        int pi = static_cast<int>((phi + M_PI) / (2.0 * M_PI) * n_lon_);
        ti = std::clamp(ti, 0, n_lat_ - 1);
        pi = ((pi % n_lon_) + n_lon_) % n_lon_;
        return static_cast<size_t>(ti) * n_lon_ + pi;
    }

    // Fixed rotation placing the galactic plane at a pleasant angle to the
    // black hole's equator (the two are unrelated in nature).
    static void galactic_to_world(const double g[3], double w[3]) {
        constexpr double ca = 0.8660254037844387, sa = 0.5;      // 30 deg tilt
        w[0] = g[0];
        w[1] = ca * g[1] - sa * g[2];
        w[2] = sa * g[1] + ca * g[2];
    }
    static void world_to_galactic(const double w[3], double g[3]) {
        constexpr double ca = 0.8660254037844387, sa = 0.5;
        g[0] = w[0];
        g[1] =  ca * w[1] + sa * w[2];
        g[2] = -sa * w[1] + ca * w[2];
    }

    // Diffuse unresolved starlight of the Galactic disc.  Peak surface
    // brightness ~21.5 mag/arcsec^2, which is what the real Milky Way has -
    // roughly 10^-12 of the solar photosphere, so it only shows up at long
    // exposure.  That contrast is a physical result, not a rendering choice.
    static spec::XYZ milky_way_radiance(const double d[3]) {
        double g[3];
        world_to_galactic(d, g);
        const double sin_b = std::clamp(g[2], -1.0, 1.0);
        const double band  = std::exp(-std::fabs(sin_b) / 0.055);
        // Longitude dependence: brighter towards the Galactic centre.
        const double l = std::atan2(g[1], g[0]);
        const double bulge = 0.55 + 0.45 * std::exp(-(l * l) / (2.0 * 0.45 * 0.45));
        // Dust lanes.
        const double lane = 0.65 + 0.35 * std::sin(l * 7.3 + sin_b * 41.0) *
                                          std::sin(l * 3.1 - sin_b * 17.0);

        constexpr double mag_per_arcsec2 = 21.5;
        constexpr double F0 = 3.24e-10;                       // W m^-2 for V = 0
        constexpr double arcsec2_per_sr = 4.2545170296152206e10;
        const double L = F0 * std::pow(10.0, -0.4 * mag_per_arcsec2) * arcsec2_per_sr;

        const double amp = L * band * bulge * std::max(0.0, lane);
        const spec::XYZ c = spec::blackbody_xyz(4800.0);      // integrated disc light
        const double inv = (c.y > 0.0) ? 1.0 / c.y : 0.0;
        return {amp * c.x * inv, amp, amp * c.z * inv};
    }
};

// ---------------------------------------------------------------------------
// The Sun, as a real body: a sphere of radius R_sun radiating as a 5772 K
// blackbody with Eddington limb darkening.  Its light is deflected by the
// black hole exactly like everything else, because it is intersected against
// the traced geodesic rather than against a straight line.
// ---------------------------------------------------------------------------
struct SunBody {
    double centre[3] = {0.0, 0.0, 0.0};   // pseudo-Cartesian, gravitational radii
    double radius    = 1.0;               // gravitational radii
    double T_eff     = phys::T_sun;

    // Gravitational redshift at the Sun's own surface: z = G M / (R c^2).
    // Tiny (2.1e-6) but free to include, and it is genuinely there.
    static constexpr double self_redshift =
        phys::G * phys::M_sun / (phys::R_sun * phys::c * phys::c);

    // Intersect the segment p0 -> p1 with the sphere.  Returns the fractional
    // position of the first hit in [0,1], or -1.
    double intersect(const double p0[3], const double p1[3]) const {
        const double d[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
        const double o[3] = {p0[0] - centre[0], p0[1] - centre[1], p0[2] - centre[2]};
        const double A = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        if (A <= 0.0) return -1.0;
        const double B = 2.0 * (o[0] * d[0] + o[1] * d[1] + o[2] * d[2]);
        const double C = o[0] * o[0] + o[1] * o[1] + o[2] * o[2] - radius * radius;
        const double disc = B * B - 4.0 * A * C;
        if (disc < 0.0) return -1.0;
        const double sq = std::sqrt(disc);
        double t0 = (-B - sq) / (2.0 * A);
        double t1 = (-B + sq) / (2.0 * A);
        if (t0 > t1) std::swap(t0, t1);
        if (t0 >= 0.0 && t0 <= 1.0) return t0;
        if (t1 >= 0.0 && t1 <= 1.0) return t1;
        return -1.0;
    }

    // Eddington's grey limb-darkening law I(mu)/I(1) = (2 + 3 mu)/5, where mu
    // is the cosine of the angle between the line of sight and the local
    // surface normal.  This is why the solar disc has a bright centre and a
    // dim rim, and it is a real radiative-transfer result, not shading.
    double limb_darkening(const double hit[3], const double view_dir[3]) const {
        double n[3] = {hit[0] - centre[0], hit[1] - centre[1], hit[2] - centre[2]};
        const double len = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        if (len <= 0.0) return 1.0;
        n[0] /= len; n[1] /= len; n[2] /= len;
        // view_dir points along the traced (past-directed) ray, i.e. into the
        // surface, so the outward cosine needs a sign flip.
        const double mu = std::clamp(-(n[0] * view_dir[0] + n[1] * view_dir[1] +
                                       n[2] * view_dir[2]), 0.0, 1.0);
        return (2.0 + 3.0 * mu) / 5.0;
    }
};

// ---------------------------------------------------------------------------
// An orbiting hot spot: a compact brightness enhancement carried around the
// disc on a circular Keplerian geodesic.
//
// This exists because a Novikov-Thorne disc is *stationary and axisymmetric*,
// which has a consequence worth being explicit about: moving the camera in
// azimuth, or letting the disc "rotate", changes nothing at all in the image.
// A film of such a scene is 240 identical frames.  Something has to genuinely
// break the symmetry, and an orbiting over-density is the standard one -
// GRAVITY has watched exactly this near Sgr A*, flares tracing loops on the
// sky with 30-70 minute periods.
//
// The spot is modelled as a Gaussian enhancement of the local emitted flux, so
// the disc stays optically thick and the "first crossing wins" rule is intact.
//
// The important part is the timing.  Each traced ray already carries the
// coordinate time it took to get here (negative, since it is traced into the
// past), so the emission time is simply
//
//     t_emit = t_observer + y[Y_T]
//
// and the spot is placed where it was *then*, not where it is now.  Light
// bending means a photon that loops around the hole arrives much later than
// one that came straight, so the secondary image shows the spot at an earlier
// orbital phase than the primary.  That echo is a real, measurable effect and
// falls straight out of integrating t along with everything else.
// ---------------------------------------------------------------------------
struct HotSpot {
    bool   enabled   = false;
    double r         = 8.0;    // orbital radius, gravitational radii
    double phi0      = 0.0;    // azimuth at t = 0
    double sigma     = 0.6;    // Gaussian radius, gravitational radii
    double contrast  = 25.0;   // peak flux enhancement over the quiescent disc
    double omega     = 0.0;    // orbital angular velocity, set by set_spin()
    int    sense     = +1;

    void set_spin(double a, int prograde) {
        sense = prograde;
        omega = prograde / (r * std::sqrt(r) + prograde * a);
    }

    // Orbital period in units of GM/c^3.
    double period() const { return (omega != 0.0) ? 2.0 * M_PI / std::fabs(omega) : 0.0; }

    // Flux enhancement factor at an equatorial point (r_hit, phi_hit) that
    // emitted at coordinate time t_emit.
    double enhancement(double r_hit, double phi_hit, double t_emit) const {
        if (!enabled) return 1.0;
        const double phi_s = phi0 + omega * t_emit;
        // Separation in the equatorial plane, measured with flat distances -
        // the spot is small compared with the radius of curvature here.
        const double dx = r_hit * std::cos(phi_hit) - r * std::cos(phi_s);
        const double dy = r_hit * std::sin(phi_hit) - r * std::sin(phi_s);
        const double d2 = dx * dx + dy * dy;
        return 1.0 + contrast * std::exp(-0.5 * d2 / (sigma * sigma));
    }
};

// ---------------------------------------------------------------------------
// Camera.  The image plane is defined in the observer's own orthonormal frame,
// so relativistic aberration and the gravitational distortion of the local sky
// are automatic consequences of using that tetrad - there is no separate
// "aberration step".
// ---------------------------------------------------------------------------
struct Camera {
    double r = 60.0;                 // Boyer-Lindquist position
    double theta = M_PI / 2.0;
    double phi = 0.0;
    double fov = 0.6;                // horizontal field of view, radians
    double yaw = 0.0;                // rotate towards +e_phi
    double pitch = 0.0;              // rotate towards -e_theta (up)
    double roll = 0.0;

    // Local orthonormal direction for a normalised image-plane coordinate
    // (sx, sy) in [-1, 1], with aspect = width/height.
    // Basis in the ZAMO frame: forward = -e_r (towards the hole),
    // right = +e_phi, up = -e_theta.
    void ray_direction(double sx, double sy, double aspect, double out[3]) const {
        const double tanf = std::tan(0.5 * fov);
        double x = sx * tanf;                 // along "right"
        double y = sy * tanf / aspect;        // along "up"
        double z = 1.0;                       // along "forward"

        // Roll about the view axis.
        if (roll != 0.0) {
            const double cr = std::cos(roll), sr = std::sin(roll);
            const double xr = cr * x - sr * y, yr = sr * x + cr * y;
            x = xr; y = yr;
        }
        // Pitch (about "right"), then yaw (about "up").
        const double cp = std::cos(pitch), sp = std::sin(pitch);
        double y2 = cp * y - sp * z;
        double z2 = sp * y + cp * z;
        const double cy = std::cos(yaw), sy2 = std::sin(yaw);
        double x2 = cy * x + sy2 * z2;
        double z3 = -sy2 * x + cy * z2;

        // Convert (right, up, forward) to the tetrad basis (e_r, e_th, e_phi).
        double v[3];
        v[0] = -z3;      // forward = -e_r
        v[1] = -y2;      // up      = -e_theta
        v[2] =  x2;      // right   = +e_phi

        const double n = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        out[0] = v[0] / n; out[1] = v[1] / n; out[2] = v[2] / n;
    }
};

} // namespace bh
