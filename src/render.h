// render.h - Backwards general-relativistic ray tracing.
#pragma once

#include "disc.h"
#include "scene.h"
#include "image.h"

#include <memory>
#include <string>

namespace bh {

struct RenderConfig {
    // Black hole
    double M_kg = 10.0 * phys::M_sun;
    double spin = 0.9;

    // Camera
    Camera cam;
    int width = 1280, height = 720;
    int sqrt_spp = 2;                  // sqrt of samples per pixel

    // Content
    bool  disc_enabled = true;
    double disc_r_out = 24.0;          // gravitational radii
    double eddington  = 0.1;           // accretion rate / Eddington rate
    int   disc_sense  = +1;
    bool  stars_enabled = true;
    int   star_count = 250000;
    uint64_t seed = 20240921ull;
    bool  sun_enabled = false;
    SunBody sun;

    // Integration
    double rtol = 1e-9;
    double atol = 1e-11;
    double r_escape = 4000.0;
    int    max_steps = 20000;

    // Execution
    int threads = 0;                   // 0 = hardware concurrency
    bool quiet = false;
};

struct RenderStats {
    long long rays = 0;
    long long steps = 0;
    long long captured = 0;
    long long disc_hits = 0;
    long long sun_hits = 0;
    long long escaped = 0;
    long long exhausted = 0;
    double max_norm_error = 0.0;       // worst |g^{ab} p_a p_b| along any ray
    double max_carter_drift = 0.0;     // worst relative drift of Carter's Q
    double seconds = 0.0;
};

// Trace the whole frame.  Returns absolute radiance in CIE XYZ.
img::Image render(const RenderConfig& cfg, RenderStats& stats);

} // namespace bh
