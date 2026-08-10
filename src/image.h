// image.h - Framebuffer, display transform, and dependency-free PNG output.
//
// The renderer produces absolute spectral radiance.  Turning that into pixels
// is a *display* problem, not a physics problem, and the two are kept strictly
// apart: nothing in here feeds back into the simulation.
#pragma once

#include "spectrum.h"

#include <string>
#include <vector>
#include <cstdint>

namespace img {

enum class ToneMap { Linear, Reinhard, ACES, Log };

struct Image {
    int width = 0, height = 0;
    std::vector<spec::XYZ> pixels;   // absolute radiance, CIE XYZ

    Image() = default;
    Image(int w, int h) : width(w), height(h), pixels(static_cast<size_t>(w) * h) {}

    spec::XYZ&       at(int x, int y)       { return pixels[static_cast<size_t>(y) * width + x]; }
    const spec::XYZ& at(int x, int y) const { return pixels[static_cast<size_t>(y) * width + x]; }
};

// Meter the scene: place the top of its histogram just below clipping.
// `key` is the tone-curve input the highlights are mapped to (1.0 lands at
// about 0.9 on screen through the ACES curve).
double auto_exposure(const Image& im, double key = 1.0);

// Veiling glare: light scattered inside a real optical system (and inside the
// human eye).  Physically it is a wide, low-amplitude PSF added to the sharp
// image; it is what makes bright sources bloom.  Optional.
void add_glare(Image& im, double strength, double radius_px);

// Apply exposure + tone curve + sRGB transfer function.  `log_decades` sets
// how many decades below middle grey the Log curve reaches; it is ignored by
// the other curves.
std::vector<uint8_t> develop(const Image& im, double exposure, ToneMap tm,
                             double log_decades = 6.0);

bool write_png(const std::string& path, const std::vector<uint8_t>& rgb, int w, int h);

// Animated PNG: a self-contained moving image that plays in any browser, with
// no external encoder. Frames must all be w*h*3 bytes. loops = 0 means forever.
bool write_apng(const std::string& path, const std::vector<std::vector<uint8_t>>& frames,
                int w, int h, int fps, int loops = 0);
bool write_ppm(const std::string& path, const std::vector<uint8_t>& rgb, int w, int h);

} // namespace img
