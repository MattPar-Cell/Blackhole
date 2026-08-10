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

    // Veiling glare: light scattered inside the optics (and inside the eye).
    // The real point spread function is heavy-tailed - Spencer et al. (1995)
    // measure roughly theta^-3 for the human eye - so this superposes a narrow
    // and a wide Gaussian rather than using a single one.
    //
    // A caveat worth knowing: any glare model has a finite kernel, and under a
    // very wide display stretch (--tonemap log with many decades) the edge of
    // that kernel becomes visible as a halo boundary around saturated sources.
    // That is a limit of the convolution, not of the physics.  Pass --glare 0
    // for wide-latitude renders.
    auto blur = [&](double sigma, std::vector<spec::XYZ>& dst) {
        const int r = std::max(1, static_cast<int>(std::ceil(3.0 * sigma)));
        std::vector<double> k(2 * r + 1);
        double sum = 0.0;
        for (int i = -r; i <= r; ++i) {
            k[i + r] = std::exp(-0.5 * (i * i) / (sigma * sigma));
            sum += k[i + r];
        }
        for (double& v : k) v /= sum;

        std::vector<spec::XYZ> tmp(im.pixels.size());
        for (int y = 0; y < h; ++y)
            for (int x = 0; x < w; ++x) {
                spec::XYZ a{};
                for (int i = -r; i <= r; ++i) {
                    const spec::XYZ& s = im.pixels[static_cast<size_t>(y) * w +
                                                   std::clamp(x + i, 0, w - 1)];
                    a.x += k[i + r] * s.x; a.y += k[i + r] * s.y; a.z += k[i + r] * s.z;
                }
                tmp[static_cast<size_t>(y) * w + x] = a;
            }
        dst.assign(im.pixels.size(), spec::XYZ{});
        for (int y = 0; y < h; ++y)
            for (int x = 0; x < w; ++x) {
                spec::XYZ a{};
                for (int i = -r; i <= r; ++i) {
                    const spec::XYZ& s =
                        tmp[static_cast<size_t>(std::clamp(y + i, 0, h - 1)) * w + x];
                    a.x += k[i + r] * s.x; a.y += k[i + r] * s.y; a.z += k[i + r] * s.z;
                }
                dst[static_cast<size_t>(y) * w + x] = a;
            }
    };

    std::vector<spec::XYZ> near_psf, far_psf;
    blur(radius_px, near_psf);
    blur(radius_px * 3.0, far_psf);

    const double core = 1.0 - strength;
    const double a = 0.7 * strength, b = 0.3 * strength;
    for (size_t i = 0; i < im.pixels.size(); ++i) {
        im.pixels[i].x = im.pixels[i].x * core + near_psf[i].x * a + far_psf[i].x * b;
        im.pixels[i].y = im.pixels[i].y * core + near_psf[i].y * a + far_psf[i].y * b;
        im.pixels[i].z = im.pixels[i].z * core + near_psf[i].z * a + far_psf[i].z * b;
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

// Filter the scanlines and wrap the result in a zlib stream: the payload of an
// IDAT (or, for animation, an fdAT) chunk.
std::vector<uint8_t> encode_image_data(const std::vector<uint8_t>& rgb, int w, int h) {
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
    return z;
}

std::vector<uint8_t> make_ihdr(int w, int h) {
    std::vector<uint8_t> ihdr;
    put_u32(ihdr, static_cast<uint32_t>(w));
    put_u32(ihdr, static_cast<uint32_t>(h));
    ihdr.push_back(8);   // bit depth
    ihdr.push_back(2);   // colour type: truecolour RGB
    ihdr.push_back(0); ihdr.push_back(0); ihdr.push_back(0);
    return ihdr;
}

bool write_file(const std::string& path, const std::vector<uint8_t>& bytes) {
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    const bool ok = std::fwrite(bytes.data(), 1, bytes.size(), f) == bytes.size();
    std::fclose(f);
    return ok;
}

const uint8_t kPngSig[8] = {137, 'P', 'N', 'G', 13, 10, 26, 10};

} // namespace

bool write_png(const std::string& path, const std::vector<uint8_t>& rgb, int w, int h) {
    if (static_cast<size_t>(w) * h * 3 != rgb.size()) return false;
    const std::vector<uint8_t> z = encode_image_data(rgb, w, h);

    std::vector<uint8_t> png;
    png.insert(png.end(), kPngSig, kPngSig + 8);
    put_chunk(png, "IHDR", make_ihdr(w, h));
    put_chunk(png, "IDAT", z);
    put_chunk(png, "IEND", {});
    return write_file(path, png);
}

// ---------------------------------------------------------------------------
// Animated PNG.  A superset of PNG (the APNG 1.0 specification): an acTL
// animation-control chunk, then per frame an fcTL frame-control chunk followed
// by the pixel data - IDAT for the first frame so that a non-animating decoder
// still shows something sensible, fdAT for the rest.  Every fcTL and fdAT
// carries a sequence number and they share one counter.
//
// This exists so there is a zero-dependency route to a moving picture: the
// result plays in any browser and in most image viewers, with no encoder to
// install.  For H.264 or anything you want to edit, write a PNG sequence and
// hand it to ffmpeg instead.
// ---------------------------------------------------------------------------
bool write_apng(const std::string& path, const std::vector<std::vector<uint8_t>>& frames,
                int w, int h, int fps, int loops) {
    if (frames.empty()) return false;
    for (const auto& f : frames)
        if (static_cast<size_t>(w) * h * 3 != f.size()) return false;

    std::vector<uint8_t> png;
    png.insert(png.end(), kPngSig, kPngSig + 8);
    put_chunk(png, "IHDR", make_ihdr(w, h));

    std::vector<uint8_t> actl;
    put_u32(actl, static_cast<uint32_t>(frames.size()));
    put_u32(actl, static_cast<uint32_t>(std::max(0, loops)));   // 0 = forever
    put_chunk(png, "acTL", actl);

    uint32_t seq = 0;
    for (size_t i = 0; i < frames.size(); ++i) {
        std::vector<uint8_t> fctl;
        put_u32(fctl, seq++);
        put_u32(fctl, static_cast<uint32_t>(w));
        put_u32(fctl, static_cast<uint32_t>(h));
        put_u32(fctl, 0);                       // x offset
        put_u32(fctl, 0);                       // y offset
        // Delay as an exact rational: 1/fps of a second.
        fctl.push_back(0); fctl.push_back(1);                          // delay_num = 1
        fctl.push_back(static_cast<uint8_t>((fps >> 8) & 0xFF));       // delay_den = fps
        fctl.push_back(static_cast<uint8_t>(fps & 0xFF));
        fctl.push_back(0);   // dispose_op = NONE
        fctl.push_back(0);   // blend_op   = SOURCE
        put_chunk(png, "fcTL", fctl);

        const std::vector<uint8_t> z = encode_image_data(frames[i], w, h);
        if (i == 0) {
            put_chunk(png, "IDAT", z);
        } else {
            std::vector<uint8_t> fdat;
            put_u32(fdat, seq++);
            fdat.insert(fdat.end(), z.begin(), z.end());
            put_chunk(png, "fdAT", fdat);
        }
    }

    put_chunk(png, "IEND", {});
    return write_file(path, png);
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
