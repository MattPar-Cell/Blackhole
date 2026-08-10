#include "image.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>

namespace img {

// ---------------------------------------------------------------------------
// Display transform
// ---------------------------------------------------------------------------

// Meter the scene the way a photographer would: expose for the highlights.
//
// These images have an enormous dynamic range - a bright accretion disc can
// outshine a sixth-magnitude star by twelve orders of magnitude - so anchoring
// on the *average* would blow the subject out completely.  Instead the top of
// the histogram is placed just below clipping, which keeps the disc's internal
// structure (the temperature gradient and the Doppler beaming asymmetry)
// spread across the tonal range.
double auto_exposure(const Image& im, double key) {
    std::vector<double> lum;
    lum.reserve(im.pixels.size());
    for (const auto& p : im.pixels)
        if (p.y > 0.0 && std::isfinite(p.y)) lum.push_back(p.y);
    if (lum.empty()) return 1.0;

    auto pct = [&](double q) {
        size_t k = static_cast<size_t>(q * (lum.size() - 1));
        k = std::min(k, lum.size() - 1);
        std::nth_element(lum.begin(), lum.begin() + k, lum.end());
        return lum[k];
    };

    const double hi = pct(0.999);     // top of the bright subject
    const double peak = pct(0.9999);  // rejects single hot pixels

    double e = 1.0;
    if (hi > 0.0) e = key / hi;
    // If the subject covers less than a thousandth of the frame the percentile
    // lands on empty sky, so fall back to keeping the true peak in range.
    if (peak > 0.0) e = std::min(e, 3.0 * key / peak);
    return e;
}

void add_glare(Image& im, double strength, double radius_px) {
    if (strength <= 0.0 || radius_px <= 0.0) return;
    const int w = im.width, h = im.height;

    // Separable Gaussian, three-pass box approximation (fast and smooth).
    const int r = std::max(1, static_cast<int>(radius_px));
    std::vector<spec::XYZ> tmp(im.pixels), buf(im.pixels.size());

    auto box_h = [&](std::vector<spec::XYZ>& src, std::vector<spec::XYZ>& dst) {
        for (int y = 0; y < h; ++y) {
            spec::XYZ acc{};
            const size_t row = static_cast<size_t>(y) * w;
            for (int x = -r; x <= r; ++x) {
                const int cx = std::clamp(x, 0, w - 1);
                acc.x += src[row + cx].x; acc.y += src[row + cx].y; acc.z += src[row + cx].z;
            }
            const double inv = 1.0 / (2 * r + 1);
            for (int x = 0; x < w; ++x) {
                dst[row + x] = {acc.x * inv, acc.y * inv, acc.z * inv};
                const int add = std::clamp(x + r + 1, 0, w - 1);
                const int sub = std::clamp(x - r, 0, w - 1);
                acc.x += src[row + add].x - src[row + sub].x;
                acc.y += src[row + add].y - src[row + sub].y;
                acc.z += src[row + add].z - src[row + sub].z;
            }
        }
    };
    auto box_v = [&](std::vector<spec::XYZ>& src, std::vector<spec::XYZ>& dst) {
        for (int x = 0; x < w; ++x) {
            spec::XYZ acc{};
            for (int y = -r; y <= r; ++y) {
                const size_t i = static_cast<size_t>(std::clamp(y, 0, h - 1)) * w + x;
                acc.x += src[i].x; acc.y += src[i].y; acc.z += src[i].z;
            }
            const double inv = 1.0 / (2 * r + 1);
            for (int y = 0; y < h; ++y) {
                dst[static_cast<size_t>(y) * w + x] = {acc.x * inv, acc.y * inv, acc.z * inv};
                const size_t ia = static_cast<size_t>(std::clamp(y + r + 1, 0, h - 1)) * w + x;
                const size_t is = static_cast<size_t>(std::clamp(y - r, 0, h - 1)) * w + x;
                acc.x += src[ia].x - src[is].x;
                acc.y += src[ia].y - src[is].y;
                acc.z += src[ia].z - src[is].z;
            }
        }
    };

    for (int pass = 0; pass < 3; ++pass) { box_h(tmp, buf); box_v(buf, tmp); }

    for (size_t i = 0; i < im.pixels.size(); ++i) {
        im.pixels[i].x = im.pixels[i].x * (1.0 - strength) + tmp[i].x * strength;
        im.pixels[i].y = im.pixels[i].y * (1.0 - strength) + tmp[i].y * strength;
        im.pixels[i].z = im.pixels[i].z * (1.0 - strength) + tmp[i].z * strength;
    }
}

// ACES filmic curve (Narkowicz's fit to the RRT+ODT).
static inline double aces(double x) {
    constexpr double a = 2.51, b = 0.03, c = 2.43, d = 0.59, e = 0.14;
    return std::clamp((x * (a * x + b)) / (x * (c * x + d) + e), 0.0, 1.0);
}

static inline double srgb_encode(double v) {
    v = std::clamp(v, 0.0, 1.0);
    return (v <= 0.0031308) ? 12.92 * v : 1.055 * std::pow(v, 1.0 / 2.4) - 0.055;
}

std::vector<uint8_t> develop(const Image& im, double exposure, ToneMap tm,
                             double log_decades) {
    const double log_k = std::pow(10.0, std::clamp(log_decades, 1.0, 20.0));
    std::vector<uint8_t> out(static_cast<size_t>(im.width) * im.height * 3);
    for (size_t i = 0; i < im.pixels.size(); ++i) {
        spec::XYZ c = im.pixels[i];
        if (!std::isfinite(c.x) || !std::isfinite(c.y) || !std::isfinite(c.z)) c = {0, 0, 0};
        c.x *= exposure; c.y *= exposure; c.z *= exposure;

        spec::RGB rgb = spec::gamut_clamp(spec::xyz_to_linear_srgb(c));
        rgb.r = std::max(rgb.r, 0.0);
        rgb.g = std::max(rgb.g, 0.0);
        rgb.b = std::max(rgb.b, 0.0);

        double v[3] = {rgb.r, rgb.g, rgb.b};
        switch (tm) {
            case ToneMap::Linear:
                for (double& x : v) x = std::clamp(x, 0.0, 1.0);
                break;
            case ToneMap::Reinhard:
                for (double& x : v) x = x / (1.0 + x);
                break;
            case ToneMap::ACES:
                for (double& x : v) x = aces(x);
                break;
            case ToneMap::Log: {
                // Wide-latitude curve, in the spirit of a long astrophotographic
                // exposure: many decades of radiance compressed onto the display
                // range so that a 10^12 contrast between an accretion disc and a
                // sixth-magnitude star can be seen at once.
                const double denom = std::log10(1.0 + log_k);
                for (double& x : v)
                    x = std::clamp(std::log10(1.0 + log_k * std::max(x, 0.0)) / denom, 0.0, 1.0);
                break;
            }
        }
        for (int ch = 0; ch < 3; ++ch)
            out[i * 3 + ch] = static_cast<uint8_t>(std::lround(255.0 * srgb_encode(v[ch])));
    }
    return out;
}

// ---------------------------------------------------------------------------
// PNG writer.  Self-contained: CRC-32, Adler-32, and a small DEFLATE encoder
// (LZ77 with a hash-chain match finder, fixed Huffman codes) so the output is
// a genuinely compressed, standards-compliant PNG with no external libraries.
// ---------------------------------------------------------------------------
namespace {

uint32_t crc_table[256];
bool crc_ready = false;

void init_crc() {
    for (uint32_t n = 0; n < 256; ++n) {
        uint32_t c = n;
        for (int k = 0; k < 8; ++k) c = (c & 1) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
        crc_table[n] = c;
    }
    crc_ready = true;
}

uint32_t crc32(const uint8_t* buf, size_t len, uint32_t crc = 0xFFFFFFFFu) {
    if (!crc_ready) init_crc();
    for (size_t i = 0; i < len; ++i) crc = crc_table[(crc ^ buf[i]) & 0xFF] ^ (crc >> 8);
    return crc;
}

uint32_t adler32(const uint8_t* buf, size_t len) {
    uint32_t a = 1, b = 0;
    for (size_t i = 0; i < len; ++i) {
        a = (a + buf[i]) % 65521;
        b = (b + a) % 65521;
    }
    return (b << 16) | a;
}

// --- bit-level output, DEFLATE convention -----------------------------------
struct BitWriter {
    std::vector<uint8_t> out;
    uint32_t bitbuf = 0;
    int bitcount = 0;

    // Extra bits and headers: least-significant bit first.
    void bits(uint32_t value, int n) {
        bitbuf |= (value & ((1u << n) - 1u)) << bitcount;
        bitcount += n;
        while (bitcount >= 8) {
            out.push_back(static_cast<uint8_t>(bitbuf & 0xFF));
            bitbuf >>= 8;
            bitcount -= 8;
        }
    }
    // Huffman codes: most-significant bit of the code first.
    void code(uint32_t value, int n) {
        for (int i = n - 1; i >= 0; --i) bits((value >> i) & 1u, 1);
    }
    void flush() {
        while (bitcount > 0) bits(0, 1);
    }
};

// Fixed Huffman literal/length code (RFC 1951 section 3.2.6).
inline void emit_literal(BitWriter& w, int sym) {
    if (sym <= 143)       w.code(0x30 + sym, 8);
    else if (sym <= 255)  w.code(0x190 + sym - 144, 9);
    else if (sym <= 279)  w.code(sym - 256, 7);
    else                  w.code(0xC0 + sym - 280, 8);
}

const uint16_t kLenBase[29] = {3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31,
                               35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258};
const uint8_t kLenExtra[29] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2,
                               3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0};
const uint16_t kDistBase[30] = {1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193,
                                257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145,
                                8193, 12289, 16385, 24577};
const uint8_t kDistExtra[30] = {0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6,
                                7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13};

std::vector<uint8_t> deflate_fixed(const std::vector<uint8_t>& data) {
    BitWriter w;
    w.bits(1, 1);   // BFINAL
    w.bits(1, 2);   // BTYPE = 01, fixed Huffman

    constexpr int kWindow = 32768;
    constexpr int kMinMatch = 3, kMaxMatch = 258;
    constexpr int kHashBits = 15, kHashSize = 1 << kHashBits;
    constexpr int kMaxChain = 128;

    std::vector<int> head(kHashSize, -1);
    std::vector<int> prev(data.size(), -1);

    auto hash3 = [&](size_t i) -> uint32_t {
        return ((static_cast<uint32_t>(data[i]) << 10) ^
                (static_cast<uint32_t>(data[i + 1]) << 5) ^
                 static_cast<uint32_t>(data[i + 2])) & (kHashSize - 1);
    };

    size_t i = 0;
    const size_t n = data.size();
    while (i < n) {
        int best_len = 0, best_dist = 0;
        if (i + 2 < n) {
            const uint32_t h = hash3(i);
            int cand = head[h];
            int chain = 0;
            const size_t max_here = std::min<size_t>(kMaxMatch, n - i);
            while (cand >= 0 && chain++ < kMaxChain) {
                const size_t dist = i - static_cast<size_t>(cand);
                if (dist == 0 || dist > kWindow) break;
                size_t l = 0;
                while (l < max_here && data[cand + l] == data[i + l]) ++l;
                if (static_cast<int>(l) > best_len) {
                    best_len = static_cast<int>(l);
                    best_dist = static_cast<int>(dist);
                    if (best_len >= static_cast<int>(max_here)) break;
                }
                cand = prev[cand];
            }
            // Insert current position into the chain.
            prev[i] = head[h];
            head[h] = static_cast<int>(i);
        }

        if (best_len >= kMinMatch) {
            int lc = 0;
            while (lc < 28 && kLenBase[lc + 1] <= best_len) ++lc;
            emit_literal(w, 257 + lc);
            if (kLenExtra[lc]) w.bits(static_cast<uint32_t>(best_len - kLenBase[lc]), kLenExtra[lc]);
            int dc = 0;
            while (dc < 29 && kDistBase[dc + 1] <= best_dist) ++dc;
            w.code(static_cast<uint32_t>(dc), 5);   // fixed distance codes: 5 bits
            if (kDistExtra[dc]) w.bits(static_cast<uint32_t>(best_dist - kDistBase[dc]), kDistExtra[dc]);

            // Register the skipped positions so later matches can find them.
            for (size_t k = i + 1; k < i + static_cast<size_t>(best_len) && k + 2 < n; ++k) {
                const uint32_t h2 = hash3(k);
                prev[k] = head[h2];
                head[h2] = static_cast<int>(k);
            }
            i += static_cast<size_t>(best_len);
        } else {
            emit_literal(w, data[i]);
            ++i;
        }
    }
    emit_literal(w, 256);   // end of block
    w.flush();
    return w.out;
}

void put_u32(std::vector<uint8_t>& v, uint32_t x) {
    v.push_back(static_cast<uint8_t>(x >> 24));
    v.push_back(static_cast<uint8_t>(x >> 16));
    v.push_back(static_cast<uint8_t>(x >> 8));
    v.push_back(static_cast<uint8_t>(x));
}

void put_chunk(std::vector<uint8_t>& v, const char tag[4], const std::vector<uint8_t>& payload) {
    put_u32(v, static_cast<uint32_t>(payload.size()));
    const size_t start = v.size();
    v.insert(v.end(), tag, tag + 4);
    v.insert(v.end(), payload.begin(), payload.end());
    const uint32_t c = crc32(v.data() + start, v.size() - start) ^ 0xFFFFFFFFu;
    put_u32(v, c);
}

} // namespace

bool write_png(const std::string& path, const std::vector<uint8_t>& rgb, int w, int h) {
    if (static_cast<size_t>(w) * h * 3 != rgb.size()) return false;

    // PNG scanline filtering.  Choosing per row with the standard minimum-sum-
    // of-absolute-differences heuristic gives the LZ77 stage much more to work
    // with on smooth astronomical images.
    const size_t stride = static_cast<size_t>(w) * 3;
    std::vector<uint8_t> raw;
    raw.reserve((stride + 1) * h);
    std::vector<uint8_t> line(stride), best(stride);

    auto paeth = [](int a, int b, int c) {
        const int p = a + b - c;
        const int pa = std::abs(p - a), pb = std::abs(p - b), pc = std::abs(p - c);
        if (pa <= pb && pa <= pc) return a;
        return (pb <= pc) ? b : c;
    };

    for (int y = 0; y < h; ++y) {
        const uint8_t* cur = rgb.data() + static_cast<size_t>(y) * stride;
        const uint8_t* up  = (y > 0) ? cur - stride : nullptr;
        int best_type = 0;
        long best_score = -1;
        for (int type = 0; type < 5; ++type) {
            long score = 0;
            for (size_t x = 0; x < stride; ++x) {
                const int a = (x >= 3) ? cur[x - 3] : 0;
                const int b = up ? up[x] : 0;
                const int c = (up && x >= 3) ? up[x - 3] : 0;
                int v = 0;
                switch (type) {
                    case 0: v = cur[x]; break;
                    case 1: v = cur[x] - a; break;
                    case 2: v = cur[x] - b; break;
                    case 3: v = cur[x] - ((a + b) >> 1); break;
                    case 4: v = cur[x] - paeth(a, b, c); break;
                }
                line[x] = static_cast<uint8_t>(v);
                score += std::abs(static_cast<int8_t>(line[x]));
            }
            if (best_score < 0 || score < best_score) {
                best_score = score;
                best_type = type;
                best = line;
            }
        }
        raw.push_back(static_cast<uint8_t>(best_type));
        raw.insert(raw.end(), best.begin(), best.end());
    }

    // zlib stream: 2-byte header, DEFLATE payload, Adler-32 of the raw data.
    std::vector<uint8_t> z;
    z.push_back(0x78); z.push_back(0x01);      // CM=8, CINFO=7, FCHECK ok, no dict
    const std::vector<uint8_t> comp = deflate_fixed(raw);
    z.insert(z.end(), comp.begin(), comp.end());
    put_u32(z, adler32(raw.data(), raw.size()));

    std::vector<uint8_t> png;
    const uint8_t sig[8] = {137, 'P', 'N', 'G', '\r', '\n', 26, '\n'};
    png.insert(png.end(), sig, sig + 8);

    std::vector<uint8_t> ihdr;
    put_u32(ihdr, static_cast<uint32_t>(w));
    put_u32(ihdr, static_cast<uint32_t>(h));
    ihdr.push_back(8);   // bit depth
    ihdr.push_back(2);   // colour type: truecolour RGB
    ihdr.push_back(0); ihdr.push_back(0); ihdr.push_back(0);
    put_chunk(png, "IHDR", ihdr);
    put_chunk(png, "IDAT", z);
    put_chunk(png, "IEND", {});

    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    const bool ok = std::fwrite(png.data(), 1, png.size(), f) == png.size();
    std::fclose(f);
    return ok;
}

bool write_ppm(const std::string& path, const std::vector<uint8_t>& rgb, int w, int h) {
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    std::fprintf(f, "P6\n%d %d\n255\n", w, h);
    const bool ok = std::fwrite(rgb.data(), 1, rgb.size(), f) == rgb.size();
    std::fclose(f);
    return ok;
}

} // namespace img
