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
    double sun_distance = 0.0;      // 0 = auto placement
    double inclination_deg = 85.0;
    std::string preset = "sgra";
};

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
  --out FILE           output file (.png or .ppm)

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
        else if (s == "--threads") o.cfg.threads = static_cast<int>(need_num(argc, argv, i));
        else if (s == "--quiet") o.cfg.quiet = true;
        else { std::fprintf(stderr, "error: unknown option '%s'\n", s.c_str()); usage(2); }
    }

    if (o.run_tests) return run_validation();

    o.cfg.cam.theta = std::clamp(o.inclination_deg, 0.5, 179.5) * M_PI / 180.0;
    o.ppm = o.out.size() > 4 && o.out.substr(o.out.size() - 4) == ".ppm";

    if (o.cfg.sun_enabled) place_sun(o);

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

    bh::RenderStats stats;
    img::Image im = bh::render(o.cfg, stats);

    if (o.glare > 0.0) img::add_glare(im, o.glare, o.glare_radius);

    const double exposure = (o.exposure > 0.0) ? o.exposure : img::auto_exposure(im);
    const std::vector<uint8_t> rgb = img::develop(im, exposure, o.tonemap, o.log_decades);

    const bool ok = o.ppm ? img::write_ppm(o.out, rgb, o.cfg.width, o.cfg.height)
                          : img::write_png(o.out, rgb, o.cfg.width, o.cfg.height);
    if (!ok) { std::fprintf(stderr, "error: could not write %s\n", o.out.c_str()); return 1; }

    if (!o.cfg.quiet) {
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
        std::printf("  exposure              %.4e%s\n", exposure, o.exposure > 0.0 ? "" : "  (auto)");
        std::printf("\nWrote %s\n", o.out.c_str());
    }
    return 0;
}
