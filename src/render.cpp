#include "render.h"

#include <atomic>
#include <chrono>
#include <cstdio>
#include <thread>
#include <vector>

namespace bh {

namespace {

// Coordinate velocity in the pseudo-Cartesian chart, from the BL velocities.
// Used to hand an escaping ray a direction on the sky at infinity.
void cartesian_velocity(double a, const State& y, const State& dy, double v[3]) {
    const double r = y[Y_R], th = y[Y_TH], ph = y[Y_PH];
    const double R = std::sqrt(r * r + a * a);
    const double sn = std::sin(th), cs = std::cos(th);
    const double sp = std::sin(ph), cp = std::cos(ph);
    const double dr = dy[Y_R], dth = dy[Y_TH], dph = dy[Y_PH];
    const double dR = (R > 0.0) ? r / R * dr : dr;

    v[0] = dR * sn * cp + R * cs * dth * cp - R * sn * sp * dph;
    v[1] = dR * sn * sp + R * cs * dth * sp + R * sn * cp * dph;
    v[2] = dr * cs - r * sn * dth;
}

// Redshift factor g = nu_obs / nu_emit for the integrated (past-directed)
// momentum k and an emitter with four-velocity u.  The camera normalisation
// k.u_cam = 1 is built into make_photon, so g reduces to 1 / (k.u_emit).
double redshift_factor(const State& y, const std::array<double, 4>& u) {
    const double ku = y[Y_PT] * u[0] + y[Y_PR] * u[1] + y[Y_PTH] * u[2] + y[Y_PPH] * u[3];
    return (ku > 1e-12) ? 1.0 / ku : 0.0;
}

struct Hit {
    enum class What { Horizon, Disc, Sun, Star, Sky, Exhausted } what = What::Exhausted;
    State y{};
    // Optically thin emission picked up on the way, in W m^-2 sr^-1 Hz^-1 at
    // the jet's reference frequency:  Integral j_nu0 g^(2+alpha) dlambda.
    double jet = 0.0;
};

// Four-point Gauss-Legendre nodes and weights on [0, 1].  The jet emissivity
// is smooth by construction (soft ends, Gaussian transverse profile), so four
// samples per accepted step resolve it without banding.
constexpr double kGLx[4] = {0.06943184420297371, 0.33000947820757187,
                            0.66999052179242813, 0.93056815579702629};
constexpr double kGLw[4] = {0.17392742256872693, 0.32607257743127307,
                            0.32607257743127307, 0.17392742256872693};

// Shortest distance from the origin to the segment p0 -> p1.  Used to skip the
// volumetric sampling entirely for the great majority of steps, which never go
// anywhere near the jet.
double segment_distance_to_origin(const double p0[3], const double p1[3]) {
    const double d[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
    const double dd = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    double t = 0.0;
    if (dd > 0.0) t = std::clamp(-(p0[0] * d[0] + p0[1] * d[1] + p0[2] * d[2]) / dd, 0.0, 1.0);
    const double q[3] = {p0[0] + t * d[0], p0[1] + t * d[1], p0[2] + t * d[2]};
    return std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
}

} // namespace

// ---------------------------------------------------------------------------
// Trace one photon backwards from the camera.
// ---------------------------------------------------------------------------
static Hit trace(const RenderConfig& cfg, const NovikovThorneDisc* disc,
                 const double dir[3], const Tetrad& tet, RenderStats& local) {
    const double a = cfg.spin;
    const double r_h = horizon_outer(a);
    const double r_stop = r_h * 1.0 + 1e-4 + (r_h - horizon_inner(a)) * 1e-3;

    State y = make_photon(a, cfg.cam.r, cfg.cam.theta, cfg.cam.phi, tet, dir, true);

    const Invariants iv0 = invariants(a, y);
    const double Q0 = std::fabs(iv0.Q) + 1e-12;

    Stepper st(a, cfg.rtol, cfg.atol);
    st.reset();

    // The momentum itself is past-directed (dt/dl < 0), so the affine
    // parameter still advances forwards: increasing lambda walks the photon
    // backwards in time, away from the camera and towards its source.
    double h = 0.02 * std::max(1.0, cfg.cam.r);

    Hit hit;
    hit.y = y;

    double p0[3];
    bl_to_cartesian(a, y[Y_R], y[Y_TH], y[Y_PH], p0);

    // -----------------------------------------------------------------------
    // Volumetric emission.  The jet is optically thin, so unlike the disc or a
    // stellar surface it does not stop the ray: its light is *added* to
    // whatever the ray eventually finds behind it.  The transfer integral
    //
    //     I_nu(obs) / (nu_obs/nu_0)^-alpha  =  Integral j_nu0 g^(2+alpha) dlambda
    //
    // is accumulated with four-point Gauss-Legendre quadrature on each accepted
    // step, evaluated on the dense-output polynomial so the samples lie on the
    // true geodesic rather than on a chord.  `s_end` truncates the step at an
    // occluding surface, so a jet behind the disc is correctly hidden by it.
    // -----------------------------------------------------------------------
    const Jet& jet = cfg.jet;
    auto sample_jet = [&](const DenseSegment& dense, double h_step, double s_end) {
        if (!jet.enabled || s_end <= 0.0) return;
        double sum = 0.0;
        for (int q = 0; q < 4; ++q) {
            const State ys = dense.eval(kGLx[q] * s_end);
            double z, rho;
            Jet::cylindrical(a, ys[Y_R], ys[Y_TH], z, rho);
            const double sh = jet.shape(z, rho);
            if (sh <= 0.0) continue;
            const Geom gg = geom_at(a, ys[Y_R], ys[Y_TH]);
            const double gf = redshift_factor(ys, jet.four_velocity(gg));
            if (gf <= 0.0) continue;
            sum += kGLw[q] * sh * jet.boost(gf);
        }
        hit.jet += jet.j_nu0 * sum * h_step * s_end;
    };

    for (int step = 0; step < cfg.max_steps; ++step) {
        State ynew{};
        DenseSegment dense;

        // Limit the step so the geodesic stays well resolved near the hole and
        // so that no single step can skip over the disc plane.
        const double hmax = 0.05 * std::max(1.0, y[Y_R] - r_h) + 0.01;
        if (h > hmax) h = hmax;

        int tries = 0;
        while (!st.try_step(y, h, ynew, dense)) {
            if (++tries > 40 || std::fabs(h) < 1e-14) { hit.what = Hit::What::Horizon; return hit; }
        }
        ++local.steps;

        // --- accuracy monitors ------------------------------------------------
        {
            const Invariants iv = invariants(a, ynew);
            local.max_norm_error = std::max(local.max_norm_error, std::fabs(iv.norm_rel));
            local.max_carter_drift =
                std::max(local.max_carter_drift, std::fabs(iv.Q - iv0.Q) / Q0);
        }

        double p1[3];
        bl_to_cartesian(a, ynew[Y_R], ynew[Y_TH], ynew[Y_PH], p1);

        // Most steps are nowhere near the jet; one distance test skips them.
        const bool near_jet =
            jet.enabled && segment_distance_to_origin(p0, p1) < jet.r_bound;

        // --- event: the neutron star surface ----------------------------------
        // There is no horizon to fall through: light stops here.  The surface
        // is a sphere of constant Schwarzschild r, so the crossing is found by
        // bisecting r on the dense-output polynomial.
        if (cfg.star.enabled && ynew[Y_R] <= cfg.star.R) {
            double lo = 0.0, hi = 1.0;
            for (int k = 0; k < 60; ++k) {
                const double mid = 0.5 * (lo + hi);
                if (dense.eval(mid)[Y_R] <= cfg.star.R) hi = mid; else lo = mid;
            }
            const double s_end = 0.5 * (lo + hi);
            if (near_jet) sample_jet(dense, dense.h, s_end);
            hit.what = Hit::What::Star;
            hit.y = dense.eval(s_end);
            return hit;
        }

        // --- event: swallowed by the hole -------------------------------------
        if (ynew[Y_R] <= r_stop) {
            if (near_jet) sample_jet(dense, dense.h, 1.0);
            hit.what = Hit::What::Horizon;
            hit.y = ynew;
            return hit;
        }

        // --- event: equatorial plane crossing (the disc) ----------------------
        if (disc) {
            const double c0 = std::cos(y[Y_TH]);
            const double c1 = std::cos(ynew[Y_TH]);
            if (c0 * c1 < 0.0) {
                double lo = 0.0, hi = 1.0;
                for (int k = 0; k < 60; ++k) {
                    const double mid = 0.5 * (lo + hi);
                    const double cm = std::cos(dense.eval(mid)[Y_TH]);
                    if (cm * c0 <= 0.0) hi = mid; else lo = mid;
                }
                const double s_end = 0.5 * (lo + hi);
                const State yc = dense.eval(s_end);
                if (disc->contains(yc[Y_R])) {
                    if (near_jet) sample_jet(dense, dense.h, s_end);
                    hit.what = Hit::What::Disc;
                    hit.y = yc;
                    return hit;
                }
            }
        }

        // --- event: the Sun ---------------------------------------------------
        if (cfg.sun_enabled) {
            const double* C = cfg.sun.centre;
            const double seg[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
            const double seg_len = std::sqrt(seg[0] * seg[0] + seg[1] * seg[1] + seg[2] * seg[2]);
            const double dc = std::sqrt((p0[0] - C[0]) * (p0[0] - C[0]) +
                                        (p0[1] - C[1]) * (p0[1] - C[1]) +
                                        (p0[2] - C[2]) * (p0[2] - C[2]));
            if (dc <= seg_len + cfg.sun.radius * 1.5) {
                // Sample the true geodesic across the step rather than trusting
                // the chord, then bisect on the first bracketed crossing.
                auto f = [&](double s) {
                    const State ys = dense.eval(s);
                    double q[3];
                    bl_to_cartesian(a, ys[Y_R], ys[Y_TH], ys[Y_PH], q);
                    const double dx = q[0] - C[0], dy2 = q[1] - C[1], dz = q[2] - C[2];
                    return dx * dx + dy2 * dy2 + dz * dz - cfg.sun.radius * cfg.sun.radius;
                };
                constexpr int kProbe = 12;
                if (f(0.0) <= 0.0) {   // camera already inside: treat as a hit
                    hit.what = Hit::What::Sun; hit.y = y; return hit;
                }
                double prev_s = 0.0;
                for (int k = 1; k <= kProbe; ++k) {
                    const double s = static_cast<double>(k) / kProbe;
                    if (f(s) <= 0.0) {
                        double lo = prev_s, hi2 = s;
                        for (int b = 0; b < 50; ++b) {
                            const double mid = 0.5 * (lo + hi2);
                            if (f(mid) <= 0.0) hi2 = mid; else lo = mid;
                        }
                        const double s_end = 0.5 * (lo + hi2);
                        if (near_jet) sample_jet(dense, dense.h, s_end);
                        hit.what = Hit::What::Sun;
                        hit.y = dense.eval(s_end);
                        return hit;
                    }
                    prev_s = s;
                }
            }
        }

        // Nothing stopped the ray inside this step, so the whole of it counts.
        if (near_jet) sample_jet(dense, dense.h, 1.0);

        // --- event: escape to the sky ----------------------------------------
        if (ynew[Y_R] > cfg.r_escape) {
            hit.what = Hit::What::Sky;
            hit.y = ynew;
            return hit;
        }

        y = ynew;
        p0[0] = p1[0]; p0[1] = p1[1]; p0[2] = p1[2];
    }

    hit.what = Hit::What::Exhausted;
    hit.y = y;
    return hit;
}

// ---------------------------------------------------------------------------
// Turn a hit into radiance.
// ---------------------------------------------------------------------------
static spec::XYZ shade_hit(const RenderConfig& cfg, const NovikovThorneDisc* disc,
                           const StarField* stars, const Hit& hit, double psf_sigma,
                           RenderStats& local) {
    const double a = cfg.spin;

    switch (hit.what) {
        case Hit::What::Horizon:
            ++local.captured;
            return {0.0, 0.0, 0.0};    // nothing comes back out. Ever.

        case Hit::What::Disc: {
            ++local.disc_hits;
            const double r = hit.y[Y_R];
            double T_em = disc->temperature(r);
            if (T_em <= 0.0) return {0.0, 0.0, 0.0};

            // An orbiting hot spot enhances the locally emitted flux. Since
            // F = sigma T^4, a factor f on the flux is f^{1/4} on the
            // temperature. The emission time is the observer's time plus the
            // (negative) coordinate time the ray spent getting here, so the
            // spot is sampled where it was when the light left, not where it
            // is now.
            if (cfg.hotspot.enabled) {
                const double t_emit = cfg.observer_time + hit.y[Y_T];
                const double f = cfg.hotspot.enhancement(r, hit.y[Y_PH], t_emit);
                T_em *= std::pow(f, 0.25);
            }
            const std::array<double, 4> u = disc->four_velocity(r);
            const double g = redshift_factor(hit.y, u);
            if (g <= 0.0) return {0.0, 0.0, 0.0};
            // Liouville: a blackbody stays a blackbody, at temperature g*T.
            // Gravitational redshift, transverse and longitudinal Doppler, and
            // relativistic beaming are all inside this single factor.
            return spec::blackbody_xyz(g * T_em);
        }

        case Hit::What::Sun: {
            ++local.sun_hits;
            // The Sun is at rest in Boyer-Lindquist coordinates.
            const Geom gm = geom_at(a, hit.y[Y_R], hit.y[Y_TH]);
            const MetricLower ml = metric_lower(gm);
            if (ml.tt >= 0.0) return {0.0, 0.0, 0.0};   // inside the ergosphere: no static observer
            const std::array<double, 4> u = {1.0 / std::sqrt(-ml.tt), 0.0, 0.0, 0.0};
            const double g = redshift_factor(hit.y, u);
            if (g <= 0.0) return {0.0, 0.0, 0.0};

            double p[3];
            bl_to_cartesian(a, hit.y[Y_R], hit.y[Y_TH], hit.y[Y_PH], p);
            const State d = geodesic_rhs(a, hit.y);
            double v[3];
            cartesian_velocity(a, hit.y, d, v);
            const double vn = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            if (vn > 0.0) { v[0] /= vn; v[1] /= vn; v[2] /= vn; }

            const double limb = cfg.sun.limb_darkening(p, v);
            // The photosphere's own gravitational redshift is included too.
            const double T_obs = g * cfg.sun.T_eff / (1.0 + SunBody::self_redshift);
            spec::XYZ c = spec::blackbody_xyz(T_obs);
            c.x *= limb; c.y *= limb; c.z *= limb;
            return c;
        }

        case Hit::What::Star: {
            ++local.star_hits;
            const NeutronStar& st = cfg.star;
            const double r = hit.y[Y_R], th = hit.y[Y_TH], ph = hit.y[Y_PH];

            // The surface rotates rigidly; its four-velocity carries the
            // Doppler shift, aberration and beaming all at once.
            const std::array<double, 4> u = st.surface_four_velocity(th);
            const double ku = hit.y[Y_PT] * u[0] + hit.y[Y_PPH] * u[3];
            if (!(ku > 1e-12)) return {0.0, 0.0, 0.0};
            const double g = 1.0 / ku;

            // Where the caps were when the light left, not where they are now.
            const double t_emit = cfg.observer_time + hit.y[Y_T];
            const double T_local = st.temperature_at(th, ph, t_emit);

            // Angle to the local normal, in the frame of the moving surface.
            // The radial leg of the tetrad is unaffected by a purely azimuthal
            // boost, so it is the static one: e_(r)^mu = (0, sqrt(1-2M/r),0,0).
            const double f = 1.0 - 2.0 / r;
            const double mu = (f > 0.0)
                ? std::clamp(-hit.y[Y_PR] * std::sqrt(f) / ku, 0.0, 1.0)
                : 1.0;
            // Eddington grey atmosphere, the same law as the solar photosphere.
            const double limb = (2.0 + 3.0 * mu) / 5.0;

            spec::XYZ c = spec::blackbody_xyz(g * T_local);
            c.x *= limb; c.y *= limb; c.z *= limb;
            return c;
        }

        case Hit::What::Sky: {
            ++local.escaped;
            if (!stars) return {0.0, 0.0, 0.0};
            const State d = geodesic_rhs(a, hit.y);
            double v[3];
            cartesian_velocity(a, hit.y, d, v);
            const double vn = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            if (vn <= 0.0) return {0.0, 0.0, 0.0};
            v[0] /= vn; v[1] /= vn; v[2] /= vn;
            return stars->radiance(v, psf_sigma);
        }

        case Hit::What::Exhausted:
        default:
            ++local.exhausted;
            return {0.0, 0.0, 0.0};
    }
}

// Whatever the ray ended on, add the optically thin jet emission it collected
// on the way there.  `jet_colour` is the CIE colour of a unit-amplitude
// synchrotron power law: the spectrum is scale-free under Doppler and
// gravitational shifts, so the whole jet is one colour and only the amplitude
// carries the relativistic physics.
static spec::XYZ shade(const RenderConfig& cfg, const NovikovThorneDisc* disc,
                       const StarField* stars, const Hit& hit, double psf_sigma,
                       const spec::XYZ& jet_colour, RenderStats& local) {
    spec::XYZ c = shade_hit(cfg, disc, stars, hit, psf_sigma, local);
    if (hit.jet > 0.0) {
        ++local.jet_rays;
        c.x += jet_colour.x * hit.jet;
        c.y += jet_colour.y * hit.jet;
        c.z += jet_colour.z * hit.jet;
    }
    return c;
}

// ---------------------------------------------------------------------------
img::Image render(const RenderConfig& cfg, RenderStats& stats) {
    const auto t_start = std::chrono::steady_clock::now();

    img::Image image(cfg.width, cfg.height);

    std::unique_ptr<NovikovThorneDisc> disc;
    if (cfg.disc_enabled)
        disc = std::make_unique<NovikovThorneDisc>(cfg.M_kg, cfg.spin, cfg.disc_r_out,
                                                   cfg.eddington, cfg.disc_sense);
    std::unique_ptr<StarField> stars;
    if (cfg.stars_enabled)
        stars = std::make_unique<StarField>(cfg.star_count, cfg.seed);

    const Geom gcam = geom_at(cfg.spin, cfg.cam.r, cfg.cam.theta);
    const Tetrad tet = zamo_tetrad(gcam);

    const double aspect = static_cast<double>(cfg.width) / cfg.height;
    // Angular size of one pixel, used as the star point-spread width.
    const double px_ang = cfg.cam.fov / cfg.width;
    const double psf_sigma = 0.55 * px_ang;

    const spec::XYZ jet_colour = cfg.jet.enabled ? spec::powerlaw_xyz(cfg.jet.alpha)
                                                 : spec::XYZ{};

    int nthreads = cfg.threads > 0 ? cfg.threads
                                   : static_cast<int>(std::thread::hardware_concurrency());
    if (nthreads <= 0) nthreads = 1;

    std::atomic<int> next_row{0};
    std::atomic<int> rows_done{0};
    std::vector<RenderStats> per_thread(nthreads);

    auto worker = [&](int tid) {
        RenderStats& local = per_thread[tid];
        const int S = std::max(1, cfg.sqrt_spp);
        const double inv_S = 1.0 / S;

        for (;;) {
            const int y = next_row.fetch_add(1);
            if (y >= cfg.height) break;

            for (int x = 0; x < cfg.width; ++x) {
                // Deterministic per-pixel jitter: identical output regardless
                // of thread count or scheduling.
                Rng rng(cfg.seed ^ (static_cast<uint64_t>(y) << 32) ^
                        (static_cast<uint64_t>(x) * 0x9E3779B97F4A7C15ull));
                spec::XYZ acc{};
                for (int sy = 0; sy < S; ++sy) {
                    for (int sx = 0; sx < S; ++sx) {
                        const double jx = (sx + rng.uniform()) * inv_S;
                        const double jy = (sy + rng.uniform()) * inv_S;
                        const double u = (2.0 * (x + jx) / cfg.width - 1.0);
                        const double v = (1.0 - 2.0 * (y + jy) / cfg.height);

                        double dir[3];
                        cfg.cam.ray_direction(u, v, aspect, dir);

                        ++local.rays;
                        const Hit hit = trace(cfg, disc.get(), dir, tet, local);
                        const spec::XYZ c = shade(cfg, disc.get(), stars.get(), hit,
                                                  psf_sigma, jet_colour, local);
                        acc.x += c.x; acc.y += c.y; acc.z += c.z;
                    }
                }
                const double w = 1.0 / (S * S);
                image.at(x, y) = {acc.x * w, acc.y * w, acc.z * w};
                local.total_flux += acc.y * w;
            }

            const int done = rows_done.fetch_add(1) + 1;
            if (!cfg.quiet && tid == 0 && (done % 8 == 0 || done == cfg.height)) {
                std::fprintf(stderr, "\r  tracing %5.1f%%  (%d/%d rows)",
                             100.0 * done / cfg.height, done, cfg.height);
                std::fflush(stderr);
            }
        }
    };

    std::vector<std::thread> pool;
    pool.reserve(nthreads);
    for (int t = 0; t < nthreads; ++t) pool.emplace_back(worker, t);
    for (auto& th : pool) th.join();
    if (!cfg.quiet) std::fprintf(stderr, "\r  tracing 100.0%%  (%d/%d rows)\n",
                                 cfg.height, cfg.height);

    for (const auto& s : per_thread) {
        stats.rays += s.rays;
        stats.steps += s.steps;
        stats.captured += s.captured;
        stats.disc_hits += s.disc_hits;
        stats.sun_hits += s.sun_hits;
        stats.star_hits += s.star_hits;
        stats.jet_rays += s.jet_rays;
        stats.total_flux += s.total_flux;
        stats.escaped += s.escaped;
        stats.exhausted += s.exhausted;
        stats.max_norm_error = std::max(stats.max_norm_error, s.max_norm_error);
        stats.max_carter_drift = std::max(stats.max_carter_drift, s.max_carter_drift);
    }
    stats.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
    return image;
}

} // namespace bh
