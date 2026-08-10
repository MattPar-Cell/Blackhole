// main.cpp - command line driver.

#include "render.h"
#include "constants.h"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include <thread>
#include <algorithm>
#include <chrono>

int run_validation();   // validate.cpp

namespace {

struct Options {
    bh::RenderConfig cfg;
    std::string out = "blackhole.png";
    img::ToneMap tonemap = img::ToneMap::ACES;
    double exposure = 0.0;          // 0 = auto
    double glare = 0.12;
    double log_decades = 6.0;
    double glare_radius = 3.0;
    bool run_tests = false;
    bool ppm = false;

    bool user_set_fov = false;
    bool user_set_distance = false;

    // Animation
    int    frames = 0;              // 0 = single still image
    int    fps = 24;
    double duration = 0.0;          // seconds of video; overrides `frames`
    double incl_to = -1e9;          // end-of-sequence inclination, degrees
    double dist_to = -1.0;          // end-of-sequence camera radius
    double fov_to  = -1.0;          // end-of-sequence field of view, radians
    double orbits  = 1.0;           // hot-spot orbits covered by the sequence
    bool   apng = false;
    std::string frame_pattern;      // printf pattern for a PNG sequence
    double time_span = 0.0;         // coordinate time covered, in GM/c^3

    double sun_distance = 0.0;      // 0 = auto placement
    double inclination_deg = 85.0;
    std::string preset = "sgra";
};

// Interpolate the camera and the clock to fraction u in [0, 1] of the sequence.
void set_frame(const Options& o, bh::RenderConfig& fc, double u) {
    if (o.incl_to > -1e8)
        fc.cam.theta = std::clamp(o.inclination_deg + u * (o.incl_to - o.inclination_deg),
                                  0.5, 179.5) * M_PI / 180.0;
    if (o.dist_to > 0.0) fc.cam.r = o.cfg.cam.r + u * (o.dist_to - o.cfg.cam.r);
    if (o.fov_to  > 0.0) fc.cam.fov = o.cfg.cam.fov + u * (o.fov_to - o.cfg.cam.fov);
    fc.observer_time = u * o.time_span;
}

void accumulate(bh::RenderStats& a, const bh::RenderStats& b) {
    a.rays += b.rays; a.steps += b.steps; a.captured += b.captured;
    a.disc_hits += b.disc_hits; a.sun_hits += b.sun_hits;
    a.escaped += b.escaped; a.exhausted += b.exhausted;
    a.seconds += b.seconds;
    a.max_norm_error = std::max(a.max_norm_error, b.max_norm_error);
    a.max_carter_drift = std::max(a.max_carter_drift, b.max_carter_drift);
}

void print_stats(const Options& o, const bh::RenderStats& stats, double exposure) {
    std::printf("Rendered in %.2f s  (%.3g rays, %.3g integration steps, %.1f steps/ray)\n",
                stats.seconds, static_cast<double>(stats.rays),
                static_cast<double>(stats.steps),
                stats.rays ? static_cast<double>(stats.steps) / stats.rays : 0.0);
    std::printf("  captured by the hole  %5.1f%%\n", 100.0 * stats.captured / std::max(1LL, stats.rays));
    std::printf("  hit the disc          %5.1f%%\n", 100.0 * stats.disc_hits / std::max(1LL, stats.rays));
    if (o.cfg.sun_enabled)
        std::printf("  hit the Sun           %5.1f%%\n", 100.0 * stats.sun_hits / std::max(1LL, stats.rays));
    std::printf("  escaped to the sky    %5.1f%%\n", 100.0 * stats.escaped / std::max(1LL, stats.rays));
    if (stats.exhausted)
        std::printf("  step limit reached    %5.1f%%\n", 100.0 * stats.exhausted / std::max(1LL, stats.rays));
    std::printf("  worst |g^ab p_a p_b|  %.2e   (photons should stay exactly null)\n",
                stats.max_norm_error);
    std::printf("  worst drift in Q      %.2e   (Carter's constant)\n", stats.max_carter_drift);
    std::printf("  exposure              %.4e%s\n", exposure, o.exposure > 0.0 ? "" : "  (auto, locked)");
}

[[noreturn]] void usage(int code) {
    std::printf(R"(
blackhole - a general-relativistic ray tracer for the Kerr spacetime

USAGE
  blackhole [options]

SCENE
  --preset NAME        sgra | m87 | stellar | gargantua | custom   (default: sgra)
  --mass MSUN          black hole mass in solar masses
  --spin A             dimensionless spin a/M in [-0.9999, 0.9999]
  --distance R         camera distance in gravitational radii GM/c^2
  --inclination DEG    viewing angle from the spin axis (90 = edge on)
  --fov DEG            horizontal field of view

DISC
  --no-disc            switch the accretion disc off
  --disc-outer R       outer radius in gravitational radii
  --eddington F        accretion rate as a fraction of the Eddington rate
  --retrograde         disc counter-rotates with respect to the hole

THE SUN
  --sun                add the Sun to the frame, at the same distance as the
                       black hole, so the two angular sizes are directly and
                       honestly comparable
  --sun-distance R     place the Sun this far from the hole instead (in GM/c^2)

SKY
  --no-stars           empty background instead of a lensed star field
  --stars N            number of stars (default 250000)
  --seed N             random seed for the star field

IMAGE
  --width N  --height N
  --spp N              samples per pixel is N*N (default 2 -> 4 samples)
  --exposure X         exposure multiplier; omit for automatic
  --tonemap NAME       aces | reinhard | log | linear
  --glare X            veiling-glare strength in [0,1] (default 0.12)
  --log-decades N      decades of radiance the "log" curve spans (default 6)
  --out FILE           output file (.png, .ppm, or .apng for animation)

ANIMATION
  --duration S         make a film S seconds long (with --fps sets the count)
  --frames N           or give the frame count directly
  --fps N              frames per second (default 24)
  --hotspot R          orbiting hot spot at radius R: the thing that actually
                       moves.  A Novikov-Thorne disc is stationary and
                       axisymmetric, so without this every frame is identical
  --hotspot-contrast X peak brightness over the quiescent disc (default 25)
  --hotspot-size S     Gaussian radius in GM/c^2 (default 0.6)
  --orbits N           hot-spot orbits covered by the sequence (default 1)
  --inclination-to D   sweep the viewing angle to D degrees
  --distance-to R      sweep the camera distance to R
  --fov-to D           sweep the field of view to D degrees
                       Writing to a .apng gives a self-contained movie; any
                       other extension writes a numbered frame sequence and
                       prints the ffmpeg command to turn it into an mp4.

OTHER
  --threads N          worker threads (default: all cores)
  --test               run the physics validation suite and exit
  --quiet
  --help

EXAMPLES
  blackhole --test
  blackhole --preset sgra --sun --width 1600 --height 900 --spp 3
  blackhole --preset stellar --inclination 80 --spin 0.99
  blackhole --preset m87 --inclination 17 --tonemap log
)");
    std::exit(code);
}

double need_num(int argc, char** argv, int& i) {
    if (i + 1 >= argc) { std::fprintf(stderr, "error: %s needs a value\n", argv[i]); std::exit(2); }
    return std::atof(argv[++i]);
}
std::string need_str(int argc, char** argv, int& i) {
    if (i + 1 >= argc) { std::fprintf(stderr, "error: %s needs a value\n", argv[i]); std::exit(2); }
    return argv[++i];
}

void apply_preset(Options& o, const std::string& name) {
    if (name == "sgra") {
        // Sagittarius A*, the black hole at the centre of the Milky Way.
        o.cfg.M_kg = 4.297e6 * phys::M_sun;
        o.cfg.spin = 0.90;
        o.cfg.cam.r = 150.0;
        o.cfg.disc_r_out = 26.0;
        o.cfg.eddington = 1e-6;          // Sgr A* is extraordinarily quiescent
        o.inclination_deg = 82.0;
        o.cfg.cam.fov = 0.55;
    } else if (name == "m87") {
        o.cfg.M_kg = 6.5e9 * phys::M_sun;
        o.cfg.spin = 0.94;
        o.cfg.cam.r = 120.0;
        o.cfg.disc_r_out = 30.0;
        o.cfg.eddington = 1e-5;
        o.inclination_deg = 17.0;        // M87's jet points nearly at us
        o.cfg.cam.fov = 0.6;
    } else if (name == "stellar") {
        // A typical X-ray binary: a 10 solar-mass hole accreting from a
        // companion star.
        o.cfg.M_kg = 10.0 * phys::M_sun;
        o.cfg.spin = 0.95;
        o.cfg.cam.r = 220.0;
        o.cfg.disc_r_out = 45.0;
        o.cfg.eddington = 0.1;
        o.inclination_deg = 78.0;
        o.cfg.cam.fov = 0.55;
    } else if (name == "gargantua") {
        // A slowly spun-up supermassive hole with a bright, cool disc: the
        // configuration that makes the lensed disc most spectacular.
        o.cfg.M_kg = 1.0e8 * phys::M_sun;
        o.cfg.spin = 0.999;
        o.cfg.cam.r = 110.0;
        o.cfg.disc_r_out = 24.0;
        o.cfg.eddington = 0.02;
        o.inclination_deg = 87.0;
        o.cfg.cam.fov = 0.62;
    } else if (name != "custom") {
        std::fprintf(stderr, "error: unknown preset '%s'\n", name.c_str());
        std::exit(2);
    }
}

// ---------------------------------------------------------------------------
// Place the Sun so that it and the black hole are the same distance from the
// camera.  That is the only placement for which comparing their apparent sizes
// in the image is a comparison of their real sizes.
// ---------------------------------------------------------------------------
void place_sun(Options& o) {
    const double r_g = phys::r_g_metres(o.cfg.M_kg);
    const double R_sun_geo = phys::R_sun / r_g;             // solar radius in GM/c^2
    const double R_shadow  = 3.0 * std::sqrt(3.0);          // shadow radius, ~5.196 M

    o.cfg.sun.radius = R_sun_geo;
    o.cfg.sun.T_eff = phys::T_sun;

    // Keep the Sun clear of both the shadow and the accretion disc.
    double sep = 1.8 * (R_sun_geo + R_shadow) + 3.0;
    if (o.cfg.disc_enabled) sep = std::max(sep, o.cfg.disc_r_out * 1.5 + R_sun_geo + 3.0);
    if (o.sun_distance > 0.0) sep = o.sun_distance;

    // The visible extent of the hole is whichever is larger, its shadow or the
    // accretion disc around it.
    const double R_hole = o.cfg.disc_enabled ? std::max(R_shadow, o.cfg.disc_r_out)
                                             : R_shadow;

    // Widen the shot until both objects fit comfortably.  For a stellar-mass
    // hole this pushes the camera a long way out - which is the honest
    // answer: the Sun is about 47000 times larger than a 10 solar-mass hole.
    for (int guard = 0; guard < 400; ++guard) {
        const double d_sun = std::hypot(o.cfg.cam.r, sep);
        const double delta = std::atan2(sep, o.cfg.cam.r);
        const double a_bh  = std::atan2(R_hole, o.cfg.cam.r);
        const double a_sun = std::atan2(R_sun_geo, d_sun);
        const double half  = 1.15 * (0.5 * delta + std::max(a_bh, a_sun));
        if (half < 0.65 || o.user_set_distance) {
            if (!o.user_set_fov) o.cfg.cam.fov = 2.0 * half;
            o.cfg.cam.yaw = 0.5 * delta;
            break;
        }
        o.cfg.cam.r *= 1.25;
    }

    // The hole lies along the camera's -r_hat; displace the Sun along +phi_hat,
    // which is perpendicular to the line of sight, so that it ends up at very
    // nearly the same distance from the camera as the hole.
    const double sp = std::sin(o.cfg.cam.phi), cp = std::cos(o.cfg.cam.phi);
    o.cfg.sun.centre[0] = -sp * sep;
    o.cfg.sun.centre[1] =  cp * sep;
    o.cfg.sun.centre[2] =  0.0;
}

void print_report(const Options& o) {
    const double M = o.cfg.M_kg;
    const double r_g = phys::r_g_metres(M);
    const double r_s = 2.0 * r_g;
    const double a = o.cfg.spin;

    std::printf("\n");
    std::printf("================================================================\n");
    std::printf("  Black hole\n");
    std::printf("================================================================\n");
    std::printf("  mass                 %.4g kg  =  %.4g solar masses\n", M, M / phys::M_sun);
    std::printf("  spin a/M             %.4f   (J = %.3e kg m^2/s)\n", a,
                a * phys::G * M * M / phys::c);
    std::printf("  gravitational radius %.6e m  =  %.4g km\n", r_g, r_g / 1000.0);
    std::printf("  Schwarzschild radius %.6e m  =  %.4g km  =  %.4g R_sun\n",
                r_s, r_s / 1000.0, r_s / phys::R_sun);
    std::printf("  outer horizon r+     %.4f M   (%.4g km)\n",
                bh::horizon_outer(a), bh::horizon_outer(a) * r_g / 1000.0);
    std::printf("  ergosphere (equator) %.4f M\n", bh::ergosphere_outer(a, M_PI / 2));
    std::printf("  photon orbit         %.4f M prograde, %.4f M retrograde\n",
                bh::photon_circular_orbit(a, +1), bh::photon_circular_orbit(a, -1));
    std::printf("  ISCO                 %.4f M   (%.4g km)\n",
                bh::isco_radius(a, o.cfg.disc_sense),
                bh::isco_radius(a, o.cfg.disc_sense) * r_g / 1000.0);
    const double R_shadow = 3.0 * std::sqrt(3.0);
    std::printf("  shadow diameter      %.4f M   =  %.4g km  =  %.4g R_sun\n",
                2 * R_shadow, 2 * R_shadow * r_g / 1000.0,
                2 * R_shadow * r_g / phys::R_sun);

    if (o.cfg.disc_enabled) {
        bh::NovikovThorneDisc d(M, a, o.cfg.disc_r_out, o.cfg.eddington, o.cfg.disc_sense);
        std::printf("\n  Accretion disc (Novikov-Thorne)\n");
        std::printf("  inner / outer edge   %.3f M  ->  %.1f M\n",
                    d.inner_radius(), d.outer_radius());
        std::printf("  radiative efficiency %.4f  (%.2f%% of mc^2; fusion manages 0.7%%)\n",
                    d.efficiency(), 100.0 * d.efficiency());
        std::printf("  accretion rate       %.4e kg/s  =  %.4g solar masses per year\n",
                    d.mdot_si(), d.mdot_si() * phys::year / phys::M_sun);
        std::printf("  luminosity           %.4e W  =  %.4g L_sun  =  %.3g L_Edd\n",
                    d.luminosity(), d.luminosity() / phys::L_sun,
                    d.luminosity() / phys::L_eddington(M));
        std::printf("  peak temperature     %.4e K\n", d.peak_temperature());
        const double lam = 2.897771955e-3 / std::max(d.peak_temperature(), 1.0);
        std::printf("  peak emission at     %.4g nm  (%s)\n", lam * 1e9,
                    lam < 1e-8 ? "X-ray" : lam < 4e-7 ? "ultraviolet"
                    : lam < 7e-7 ? "visible" : "infrared");
    }

    std::printf("\n  Camera\n");
    std::printf("  distance             %.4g M  =  %.4e m  =  %.4g AU\n",
                o.cfg.cam.r, o.cfg.cam.r * r_g, o.cfg.cam.r * r_g / phys::AU);
    std::printf("  inclination          %.2f deg from the spin axis\n", o.inclination_deg);
    std::printf("  field of view        %.3f deg horizontal\n", o.cfg.cam.fov * 180.0 / M_PI);
    std::printf("  apparent shadow      %.4g deg across\n",
                2.0 * std::atan2(R_shadow, o.cfg.cam.r) * 180.0 / M_PI);

    if (o.cfg.sun_enabled) {
        const double R_sun_geo = phys::R_sun / r_g;
        const double sep = std::sqrt(o.cfg.sun.centre[0] * o.cfg.sun.centre[0] +
                                     o.cfg.sun.centre[1] * o.cfg.sun.centre[1] +
                                     o.cfg.sun.centre[2] * o.cfg.sun.centre[2]);
        const double d_sun = std::hypot(o.cfg.cam.r, sep);
        const double a_sun = 2.0 * std::atan2(R_sun_geo, d_sun);
        const double a_bh  = 2.0 * std::atan2(R_shadow, o.cfg.cam.r);

        std::printf("\n================================================================\n");
        std::printf("  Size comparison: the Sun\n");
        std::printf("================================================================\n");
        std::printf("  solar radius         %.4e m  =  %.6g GM/c^2 of this hole\n",
                    phys::R_sun, R_sun_geo);
        std::printf("  solar temperature    %.1f K   (luminosity %.4e W)\n",
                    phys::T_sun, phys::L_sun);
        std::printf("  Sun placed           %.4g M from the hole, %.4g M from the camera\n",
                    sep, d_sun);
        std::printf("\n  radius ratio         shadow / Sun  =  %.4g\n",
                    R_shadow * r_g / phys::R_sun);
        std::printf("  volume ratio         shadow / Sun  =  %.4g\n",
                    std::pow(R_shadow * r_g / phys::R_sun, 3.0));
        std::printf("  apparent diameter    black hole %.4g deg,  Sun %.4g deg\n",
                    a_bh * 180.0 / M_PI, a_sun * 180.0 / M_PI);
        std::printf("  in the image         black hole %.1f px,  Sun %.1f px across\n",
                    a_bh / (o.cfg.cam.fov / o.cfg.width),
                    a_sun / (o.cfg.cam.fov / o.cfg.width));
        if (R_shadow * r_g > phys::R_sun)
            std::printf("\n  This black hole's shadow is %.4g times wider than the Sun.\n",
                        R_shadow * r_g / phys::R_sun);
        else
            std::printf("\n  The Sun is %.4g times wider than this black hole's shadow.\n",
                        phys::R_sun / (R_shadow * r_g));
        std::printf("  Both are the same distance from the camera, so what you see\n");
        std::printf("  in the frame is their true relative size.\n");

        const double px_bh  = a_bh / (o.cfg.cam.fov / o.cfg.width);
        const double px_sun = a_sun / (o.cfg.cam.fov / o.cfg.width);
        if (px_bh < 1.5 || px_sun < 1.5) {
            const char* small = (px_bh < px_sun) ? "black hole" : "Sun";
            std::printf("\n  Note: at this resolution the %s is smaller than a pixel.\n", small);
            std::printf("  That is the honest answer, not a rendering failure - the two\n");
            std::printf("  differ in size by a factor of %.0f.  Raise --width for a bigger\n",
                        std::max(px_bh, px_sun) / std::max(std::min(px_bh, px_sun), 1e-30));
            std::printf("  frame, or try --preset sgra --sun, where the ratio is only 47.\n");
        }
    }
    std::printf("\n");
}

} // namespace

int main(int argc, char** argv) {
    Options o;
    std::vector<std::string> args(argv + 1, argv + argc);

    // First pass: the preset, so explicit flags can override it.
    for (size_t i = 0; i < args.size(); ++i)
        if (args[i] == "--preset" && i + 1 < args.size()) o.preset = args[i + 1];
    apply_preset(o, o.preset);

    for (int i = 1; i < argc; ++i) {
        const std::string s = argv[i];
        if (s == "--help" || s == "-h") usage(0);
        else if (s == "--test") o.run_tests = true;
        else if (s == "--preset") need_str(argc, argv, i);
        else if (s == "--mass") o.cfg.M_kg = need_num(argc, argv, i) * phys::M_sun;
        else if (s == "--spin") o.cfg.spin = std::clamp(need_num(argc, argv, i), -0.9999, 0.9999);
        else if (s == "--distance") { o.cfg.cam.r = need_num(argc, argv, i); o.user_set_distance = true; }
        else if (s == "--inclination") o.inclination_deg = need_num(argc, argv, i);
        else if (s == "--fov") { o.cfg.cam.fov = need_num(argc, argv, i) * M_PI / 180.0; o.user_set_fov = true; }
        else if (s == "--no-disc") o.cfg.disc_enabled = false;
        else if (s == "--disc-outer") o.cfg.disc_r_out = need_num(argc, argv, i);
        else if (s == "--eddington") o.cfg.eddington = need_num(argc, argv, i);
        else if (s == "--retrograde") o.cfg.disc_sense = -1;
        else if (s == "--sun") o.cfg.sun_enabled = true;
        else if (s == "--sun-distance") { o.cfg.sun_enabled = true; o.sun_distance = need_num(argc, argv, i); }
        else if (s == "--no-stars") o.cfg.stars_enabled = false;
        else if (s == "--stars") o.cfg.star_count = static_cast<int>(need_num(argc, argv, i));
        else if (s == "--seed") o.cfg.seed = static_cast<uint64_t>(need_num(argc, argv, i));
        else if (s == "--width") o.cfg.width = static_cast<int>(need_num(argc, argv, i));
        else if (s == "--height") o.cfg.height = static_cast<int>(need_num(argc, argv, i));
        else if (s == "--spp") o.cfg.sqrt_spp = static_cast<int>(need_num(argc, argv, i));
        else if (s == "--exposure") o.exposure = need_num(argc, argv, i);
        else if (s == "--glare") o.glare = need_num(argc, argv, i);
        else if (s == "--log-decades") o.log_decades = need_num(argc, argv, i);
        else if (s == "--tonemap") {
            const std::string t = need_str(argc, argv, i);
            if (t == "aces") o.tonemap = img::ToneMap::ACES;
            else if (t == "reinhard") o.tonemap = img::ToneMap::Reinhard;
            else if (t == "log") o.tonemap = img::ToneMap::Log;
            else if (t == "linear") o.tonemap = img::ToneMap::Linear;
            else { std::fprintf(stderr, "error: unknown tonemap '%s'\n", t.c_str()); return 2; }
        }
        else if (s == "--out") o.out = need_str(argc, argv, i);
        else if (s == "--frames") o.frames = static_cast<int>(need_num(argc, argv, i));
        else if (s == "--fps") o.fps = std::max(1, static_cast<int>(need_num(argc, argv, i)));
        else if (s == "--duration") o.duration = need_num(argc, argv, i);
        else if (s == "--orbits") o.orbits = need_num(argc, argv, i);
        else if (s == "--inclination-to") o.incl_to = need_num(argc, argv, i);
        else if (s == "--distance-to") o.dist_to = need_num(argc, argv, i);
        else if (s == "--fov-to") o.fov_to = need_num(argc, argv, i) * M_PI / 180.0;
        else if (s == "--hotspot") { o.cfg.hotspot.enabled = true; o.cfg.hotspot.r = need_num(argc, argv, i); }
        else if (s == "--hotspot-contrast") o.cfg.hotspot.contrast = need_num(argc, argv, i);
        else if (s == "--hotspot-size") o.cfg.hotspot.sigma = need_num(argc, argv, i);
        else if (s == "--threads") o.cfg.threads = static_cast<int>(need_num(argc, argv, i));
        else if (s == "--quiet") o.cfg.quiet = true;
        else { std::fprintf(stderr, "error: unknown option '%s'\n", s.c_str()); usage(2); }
    }

    if (o.run_tests) return run_validation();

    o.cfg.cam.theta = std::clamp(o.inclination_deg, 0.5, 179.5) * M_PI / 180.0;
    o.ppm = o.out.size() > 4 && o.out.substr(o.out.size() - 4) == ".ppm";

    if (o.cfg.sun_enabled) place_sun(o);

    // ---- animation setup ------------------------------------------------
    if (o.duration > 0.0) o.frames = std::max(1, static_cast<int>(std::lround(o.duration * o.fps)));
    if (o.frames > 1) {
        const size_t dot = o.out.find_last_of('.');
        const std::string ext = (dot == std::string::npos) ? "" : o.out.substr(dot);
        o.apng = (ext == ".apng");
        if (!o.apng) {
            // Anything else becomes a numbered PNG sequence.  A printf pattern
            // supplied by the user is honoured as-is.
            o.frame_pattern = (o.out.find('%') != std::string::npos)
                            ? o.out
                            : o.out.substr(0, dot == std::string::npos ? o.out.size() : dot)
                              + "_%04d.png";
        }
        o.cfg.hotspot.set_spin(o.cfg.spin, o.cfg.disc_sense);
        // The sequence covers a whole number of hot-spot orbits by default, so
        // that it loops seamlessly.  With no hot spot the scene is stationary
        // and time is irrelevant, so the span is left at zero.
        if (o.cfg.hotspot.enabled) o.time_span = o.orbits * o.cfg.hotspot.period();

        if (o.cfg.hotspot.enabled && !o.cfg.quiet) {
            const double r_g_over_c = phys::r_g_metres(o.cfg.M_kg) / phys::c;
            const double period_s = o.cfg.hotspot.period() * r_g_over_c;
            const double span_s = o.time_span * r_g_over_c;
            const double video_s = static_cast<double>(o.frames) / o.fps;
            std::printf("\n  Hot spot\n");
            std::printf("  orbital radius       %.2f M   (ISCO is at %.2f M)\n",
                        o.cfg.hotspot.r, bh::isco_radius(o.cfg.spin, o.cfg.disc_sense));
            char human[64];
            if (period_s < 90.0)          std::snprintf(human, sizeof human, "%.3g seconds", period_s);
            else if (period_s < 5400.0)   std::snprintf(human, sizeof human, "%.3g minutes", period_s / 60.0);
            else if (period_s < 1.728e5)  std::snprintf(human, sizeof human, "%.3g hours", period_s / 3600.0);
            else if (period_s < 3.156e7)  std::snprintf(human, sizeof human, "%.3g days", period_s / 86400.0);
            else                          std::snprintf(human, sizeof human, "%.3g years", period_s / phys::year);
            std::printf("  orbital period       %.4g GM/c^3  =  %s\n",
                        o.cfg.hotspot.period(), human);
            std::printf("  orbital speed        %.4f c (as measured by a local static observer)\n",
                        1.0 / std::sqrt(o.cfg.hotspot.r));
            std::printf("  sequence covers      %.3g orbits  =  %.4e s of real time\n",
                        o.orbits, span_s);
            std::printf("  played over          %.2f s at %d fps  ->  %.4g x %s\n",
                        video_s, o.fps,
                        span_s > video_s ? span_s / video_s : video_s / span_s,
                        span_s > video_s ? "faster than real time" : "slower than real time");
            std::printf("  light crossing time  %.4g GM/c^3 from the camera - comparable to\n",
                        o.cfg.cam.r);
            std::printf("                       the orbital period, so the lensed images lag\n");
            std::printf("                       the direct one by a visible fraction of a cycle\n");
        }

        const bool moving = o.cfg.hotspot.enabled || o.incl_to > -1e8 ||
                            o.dist_to > 0.0 || o.fov_to > 0.0;
        if (!moving) {
            std::fprintf(stderr,
                "error: nothing in this scene changes with time, so every frame would be\n"
                "       identical.  A Novikov-Thorne disc is stationary and axisymmetric.\n"
                "       Add --hotspot R for an orbiting bright spot, or sweep the camera\n"
                "       with --inclination-to / --distance-to / --fov-to.\n");
            return 2;
        }
    } else {
        o.cfg.hotspot.set_spin(o.cfg.spin, o.cfg.disc_sense);
    }

    // The escape sphere has to sit outside everything in the scene.
    double far = o.cfg.cam.r * 3.0 + 100.0;
    if (o.cfg.sun_enabled) {
        const double sep = std::sqrt(o.cfg.sun.centre[0] * o.cfg.sun.centre[0] +
                                     o.cfg.sun.centre[1] * o.cfg.sun.centre[1] +
                                     o.cfg.sun.centre[2] * o.cfg.sun.centre[2]);
        far = std::max(far, (sep + o.cfg.sun.radius) * 3.0);
    }
    o.cfg.r_escape = far;

    if (!o.cfg.quiet) print_report(o);

    if (!o.cfg.quiet) {
        std::fprintf(stderr, "Rendering %dx%d at %d spp on %d threads...\n",
                     o.cfg.width, o.cfg.height, o.cfg.sqrt_spp * o.cfg.sqrt_spp,
                     o.cfg.threads > 0 ? o.cfg.threads
                                       : static_cast<int>(std::thread::hardware_concurrency()));
    }

    // ---------------------------------------------------------------------
    // Still image
    // ---------------------------------------------------------------------
    if (o.frames <= 1) {
        bh::RenderStats stats;
        img::Image im = bh::render(o.cfg, stats);
        if (o.glare > 0.0) img::add_glare(im, o.glare, o.glare_radius);
        const double exposure = (o.exposure > 0.0) ? o.exposure : img::auto_exposure(im);
        const std::vector<uint8_t> rgb = img::develop(im, exposure, o.tonemap, o.log_decades);

        const bool ok = o.ppm ? img::write_ppm(o.out, rgb, o.cfg.width, o.cfg.height)
                              : img::write_png(o.out, rgb, o.cfg.width, o.cfg.height);
        if (!ok) { std::fprintf(stderr, "error: could not write %s\n", o.out.c_str()); return 1; }

        if (!o.cfg.quiet) {
            print_stats(o, stats, exposure);
            std::printf("\nWrote %s\n", o.out.c_str());
        }
        return 0;
    }

    // ---------------------------------------------------------------------
    // Animation
    // ---------------------------------------------------------------------
    //
    // Exposure has to be locked across the whole sequence.  Metering each
    // frame independently would make the film flicker as the brightest thing
    // in shot changes - and worse, it would hide the very brightness variation
    // an orbiting hot spot is there to show.  So the meter is run once, on a
    // cheap low-resolution probe of the middle frame, and then held.
    double exposure = o.exposure;
    if (exposure <= 0.0) {
        bh::RenderConfig probe = o.cfg;
        probe.width = std::max(64, o.cfg.width / 6);
        probe.height = std::max(36, o.cfg.height / 6);
        probe.sqrt_spp = 1;
        probe.star_count = std::min(o.cfg.star_count, 40000);
        probe.quiet = true;
        set_frame(o, probe, 0.5);
        bh::RenderStats pstats;
        img::Image pim = bh::render(probe, pstats);
        if (o.glare > 0.0) img::add_glare(pim, o.glare, o.glare_radius);
        exposure = img::auto_exposure(pim);
    }

    std::vector<std::vector<uint8_t>> apng_frames;
    if (o.apng) apng_frames.reserve(o.frames);

    bh::RenderStats total;
    const auto t0 = std::chrono::steady_clock::now();

    for (int f = 0; f < o.frames; ++f) {
        bh::RenderConfig fc = o.cfg;
        fc.quiet = true;
        // A pure hot-spot orbit is periodic, so the last frame must not repeat
        // the first: divide by the frame count, not by count-1.  A camera
        // sweep is not periodic and does need to reach its endpoint exactly.
        const bool sweeping = o.incl_to > -1e8 || o.dist_to > 0.0 || o.fov_to > 0.0;
        const double denom = sweeping ? std::max(1, o.frames - 1) : o.frames;
        const double u = static_cast<double>(f) / denom;
        set_frame(o, fc, u);

        bh::RenderStats fs;
        img::Image im = bh::render(fc, fs);
        if (o.glare > 0.0) img::add_glare(im, o.glare, o.glare_radius);
        std::vector<uint8_t> rgb = img::develop(im, exposure, o.tonemap, o.log_decades);

        accumulate(total, fs);

        if (o.apng) {
            apng_frames.push_back(std::move(rgb));
        } else {
            char name[1024];
            std::snprintf(name, sizeof name, o.frame_pattern.c_str(), f);
            if (!img::write_png(name, rgb, fc.width, fc.height)) {
                std::fprintf(stderr, "error: could not write %s\n", name);
                return 1;
            }
        }

        if (!o.cfg.quiet) {
            const double el = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t0).count();
            const double eta = (f > 0) ? el / (f + 1) * (o.frames - f - 1) : 0.0;
            std::fprintf(stderr, "\r  frame %4d/%d   %.0f s elapsed, ~%.0f s remaining   ",
                         f + 1, o.frames, el, eta);
            std::fflush(stderr);
        }
    }
    if (!o.cfg.quiet) std::fprintf(stderr, "\n");

    if (o.apng) {
        if (!img::write_apng(o.out, apng_frames, o.cfg.width, o.cfg.height, o.fps)) {
            std::fprintf(stderr, "error: could not write %s\n", o.out.c_str());
            return 1;
        }
    }

    if (!o.cfg.quiet) {
        print_stats(o, total, exposure);
        std::printf("\n");
        if (o.apng) {
            std::printf("Wrote %s  (%d frames at %d fps = %.1f s)\n",
                        o.out.c_str(), o.frames, o.fps,
                        static_cast<double>(o.frames) / o.fps);
            std::printf("Plays in any browser. For an mp4, or to edit it, render a PNG\n");
            std::printf("sequence instead (use --out frames/f_%%04d.png).\n");
        } else {
            std::printf("Wrote %d frames matching %s\n", o.frames, o.frame_pattern.c_str());
            std::printf("\nTurn them into a video with:\n");
            std::printf("  ffmpeg -framerate %d -i %s -c:v libx264 -pix_fmt yuv420p -crf 18 out.mp4\n",
                        o.fps, o.frame_pattern.c_str());
        }
    }
    return 0;
}
