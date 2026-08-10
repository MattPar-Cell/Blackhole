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
    enum class What { Horizon, Disc, Sun, Sky, Exhausted } what = What::Exhausted;
    State y{};
};

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

        // --- event: swallowed by the hole -------------------------------------
        if (ynew[Y_R] <= r_stop) {
            hit.what = Hit::What::Horizon;
            hit.y = ynew;
            return hit;
        }

        double p1[3];
        bl_to_cartesian(a, ynew[Y_R], ynew[Y_TH], ynew[Y_PH], p1);

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
                const State yc = dense.eval(0.5 * (lo + hi));
                if (disc->contains(yc[Y_R])) {
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
                        hit.what = Hit::What::Sun;
                        hit.y = dense.eval(0.5 * (lo + hi2));
                        return hit;
                    }
                    prev_s = s;
                }
            }
        }

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
static spec::XYZ shade(const RenderConfig& cfg, const NovikovThorneDisc* disc,
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
            const double T_em = disc->temperature(r);
            if (T_em <= 0.0) return {0.0, 0.0, 0.0};
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
                                                  psf_sigma, local);
                        acc.x += c.x; acc.y += c.y; acc.z += c.z;
                    }
                }
                const double w = 1.0 / (S * S);
                image.at(x, y) = {acc.x * w, acc.y * w, acc.z * w};
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
        stats.escaped += s.escaped;
        stats.exhausted += s.exhausted;
        stats.max_norm_error = std::max(stats.max_norm_error, s.max_norm_error);
        stats.max_carter_drift = std::max(stats.max_carter_drift, s.max_carter_drift);
    }
    stats.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
    return image;
}

} // namespace bh
