#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

static const double PI = 3.14159265358979323846;

struct V3 { double x, y, z; };
static inline V3 operator+(V3 a, V3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
static inline V3 operator-(V3 a, V3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
static inline V3 operator*(V3 a, double s) { return {a.x * s, a.y * s, a.z * s}; }
static inline double dot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static inline V3 cross(V3 a, V3 b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
static inline V3 unit(V3 a) { return a * (1.0 / std::sqrt(dot(a, a))); }

static inline double sq(double x) { return x * x; }
static inline double clamp01(double x) { return x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x); }
static inline double mix(double a, double b, double t) { return a + (b - a) * t; }
static inline double smooth(double a, double b, double x) {          
    double t = clamp01((x - a) / (b - a));
    return t * t * (3.0 - 2.0 * t);
}
static inline double smin(double a, double b, double k) {            
    if (k <= 0.0) return std::min(a, b);
    double h = std::max(k - std::fabs(a - b), 0.0) / k;
    return std::min(a, b) - h * h * k * 0.25;
}
static inline double smax(double a, double b, double k) { return -smin(-a, -b, k); }

static const double DONUT_R = 1.00, DONUT_TUBE = 0.36;


static const double CUP_X = -1.30, CUP_RAD = 1.10;              
static const double CUP_Z0 = -(DONUT_R + DONUT_TUBE);           
static const double CUP_Z1 = 1.14;                             
static const double CUP_ROUND = 0.12, WALL = 0.14, FLOOR_T = 0.18, PIT_ROUND = 0.10;
static const double CUP_ZMID = 0.5 * (CUP_Z0 + CUP_Z1), CUP_HALF = 0.5 * (CUP_Z1 - CUP_Z0);
static const double HANDLE_R = 0.55, HANDLE_TUBE = 0.15;
static const double HANDLE_ALPHA = 112.0 * PI / 180.0;         
static const double WALL_X = CUP_X + CUP_RAD;                   
static const double HANDLE_CX = (WALL_X - 0.07) - HANDLE_R * std::cos(HANDLE_ALPHA);


struct Shape {
    double cx, R, r;      
    double ca, sa;         
    double cr, sr;        
    double s;           
    double k;             
    double pitZ;          
};

static Shape shapeAt(double t) {
    double a = smooth(0.00, 0.60, t);       
    double b = smooth(0.08, 0.70, t);      
    double c = smooth(0.50, 0.92, t);        
    double alpha = mix(PI, HANDLE_ALPHA, a); 
    double overlap = 0.7 * (1.0 - a);       
    Shape S;
    S.cx = mix(0.0, HANDLE_CX, a);
    S.R = mix(DONUT_R, HANDLE_R, a);
    S.r = mix(DONUT_TUBE, HANDLE_TUBE, a);
    S.ca = std::cos(alpha);
    S.sa = std::sin(alpha);
    S.cr = std::cos(alpha - overlap);
    S.sr = std::sin(alpha - overlap);
    S.s = b;
    S.k = 0.22 * smooth(0.0, 0.10, b);
    S.pitZ = mix(CUP_Z1 + 0.5, CUP_Z0 + FLOOR_T, c);
    return S;
}


static double sdf(V3 p, const Shape& S) {
   
   
    double X = p.x - S.cx, Z = std::fabs(p.z);
    double rho = std::sqrt(X * X + Z * Z);
    double dCircle = std::sqrt(sq(rho - S.R) + p.y * p.y);                              
    double dEndH = std::sqrt(sq(S.R * S.ca - X) + sq(S.R * S.sa - Z) + p.y * p.y);      
    double dEndR = std::sqrt(sq(S.R * S.cr - X) + sq(S.R * S.sr - Z) + p.y * p.y);      
    bool inHandle = (S.ca * Z - S.sa * X) <= 0.0;                                       
    bool inRest = (S.cr * Z - S.sr * X) >= 0.0;                                          ?
    double dHandle = (inHandle ? dCircle : dEndH) - S.r;
    double dRest = (inRest ? dCircle : dEndR) - S.r;

    
    double rc = std::sqrt(sq(p.x - CUP_X) + p.y * p.y);                                 
    double qx = rc - (CUP_RAD - CUP_ROUND);
    double qz = std::fabs(p.z - CUP_ZMID) - (CUP_HALF - CUP_ROUND);
    double dCup = std::min(std::max(qx, qz), 0.0) + std::sqrt(sq(std::max(qx, 0.0)) + sq(std::max(qz, 0.0))) - CUP_ROUND;

   
    double d = smin(dHandle, mix(dRest, dCup, S.s), S.k);

  
    double u = rc - (CUP_RAD - WALL - PIT_ROUND);
    double w = (S.pitZ + PIT_ROUND) - p.z;
    double dPit = std::min(std::max(u, w), 0.0) + std::sqrt(sq(std::max(u, 0.0)) + sq(std::max(w, 0.0))) - PIT_ROUND;
    return smax(d, -dPit, 0.05);
}

struct Betti { long chi; int b0, b1, b2; };

static Betti measureTopology(const Shape& S) {
    const double X0 = -2.7, X1 = 1.6, Y0 = -1.5, Y1 = 1.5, Z0 = -1.6, Z1 = 1.6, H = 0.04;
    const int nx = static_cast<int>((X1 - X0) / H), ny = static_cast<int>((Y1 - Y0) / H), nz = static_cast<int>((Z1 - Z0) / H);
    auto id = [&](int i, int j, int k) { return (static_cast<std::size_t>(k) * ny + j) * nx + i; };

    std::vector<uint8_t> in(static_cast<std::size_t>(nx) * ny * nz);
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                in[id(i, j, k)] = sdf({X0 + (i + 0.5) * H, Y0 + (j + 0.5) * H, Z0 + (k + 0.5) * H}, S) < 0.0 ? 1 : 0;

 
    auto cell = [&](int i, int j, int k) -> int {
        return (i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz) ? in[id(i, j, k)] : 0;
    };
    long V = 0, E = 0, F = 0, C = 0;
    for (int k = 0; k <= nz; ++k)
        for (int j = 0; j <= ny; ++j)
            for (int i = 0; i <= nx; ++i) {
                int c000 = cell(i - 1, j - 1, k - 1), c100 = cell(i, j - 1, k - 1);
                int c010 = cell(i - 1, j, k - 1),     c110 = cell(i, j, k - 1);
                int c001 = cell(i - 1, j - 1, k),     c101 = cell(i, j - 1, k);
                int c011 = cell(i - 1, j, k),         c111 = cell(i, j, k);
                V += (c000 | c100 | c010 | c110 | c001 | c101 | c011 | c111);
                E += (c100 | c110 | c101 | c111);       
                E += (c010 | c110 | c011 | c111);      
                E += (c001 | c101 | c011 | c111);       
                F += (c011 | c111);                     
                F += (c101 | c111);                     
                F += (c110 | c111);                     
                C += c111;
            }


    auto components = [&](int want, bool diagonals) {
        std::vector<uint8_t> seen(in.size(), 0);
        std::vector<std::size_t> stack;
        int count = 0;
        for (std::size_t start = 0; start < in.size(); ++start) {
            if (in[start] != want || seen[start]) continue;
            ++count;
            seen[start] = 1;
            stack.push_back(start);
            while (!stack.empty()) {
                std::size_t s = stack.back();
                stack.pop_back();
                int i = static_cast<int>(s % nx), j = static_cast<int>((s / nx) % ny), k = static_cast<int>(s / (static_cast<std::size_t>(nx) * ny));
                for (int dk = -1; dk <= 1; ++dk)
                    for (int dj = -1; dj <= 1; ++dj)
                        for (int di = -1; di <= 1; ++di) {
                            int manhattan = std::abs(di) + std::abs(dj) + std::abs(dk);
                            if (manhattan == 0 || (!diagonals && manhattan != 1)) continue;
                            int a = i + di, b = j + dj, c = k + dk;
                            if (a < 0 || a >= nx || b < 0 || b >= ny || c < 0 || c >= nz) continue;
                            std::size_t n = id(a, b, c);
                            if (in[n] == want && !seen[n]) { seen[n] = 1; stack.push_back(n); }
                        }
            }
        }
        return count;
    };

    Betti r;
    r.chi = V - E + F - C;
    r.b0 = components(1, true);
    r.b2 = components(0, false) - 1;            
    r.b1 = r.b0 + r.b2 - static_cast<int>(r.chi);
    return r;
}


static const double FLOOR_Z = CUP_Z0;                       
static const V3 KEY = unit({-0.50, -0.35, 0.80});           
static const V3 FILL = unit({0.60, -0.70, 0.25});

struct Cam { V3 eye, fwd, right, up, target; double dist; };

static Cam cameraAt(double t) {
    double e = smooth(0.0, 1.0, t);
    Cam c;
    c.target = {mix(0.0, -0.88, e), 0.0, mix(0.0, -0.06, e)};   
    double az = mix(20.0, 40.0, t) * PI / 180.0, el = mix(22.0, 30.0, e) * PI / 180.0;
    c.dist = mix(8.2, 9.4, e);
    V3 dir = {std::sin(az) * std::cos(el), -std::cos(az) * std::cos(el), std::sin(el)};
    c.eye = c.target + dir * c.dist;
    c.fwd = unit(c.target - c.eye);
    c.right = unit(cross(c.fwd, {0.0, 0.0, 1.0}));
    c.up = cross(c.right, c.fwd);
    return c;
}

static V3 normalAt(V3 p, const Shape& S) {
    const double e = 0.002;
    const V3 k0{1, -1, -1}, k1{-1, -1, 1}, k2{-1, 1, -1}, k3{1, 1, 1};
    return unit(k0 * sdf(p + k0 * e, S) + k1 * sdf(p + k1 * e, S) + k2 * sdf(p + k2 * e, S) + k3 * sdf(p + k3 * e, S));
}

static double softShadow(V3 p, V3 toLight, const Shape& S) {
    double res = 1.0, t = 0.03;
    for (int i = 0; i < 28; ++i) {
        double h = sdf(p + toLight * t, S);
        res = std::min(res, 6.0 * h / t);
        if (res < 0.02) return 0.0;
        t += std::min(std::max(h, 0.03), 0.35);
        if (t > 6.0) break;
    }
    res = clamp01(res);
    return res * res * (3.0 - 2.0 * res);
}

struct Px { float L; bool hit; };                           

static Px trace(const Cam& cam, double u, double v, const Shape& S) {
    V3 rd = unit(cam.fwd + cam.right * u + cam.up * v);
    V3 ro = cam.eye;
    double tFloor = rd.z < -1e-9 ? (FLOOR_Z - ro.z) / rd.z : 1e30;

    
    const double BOUND = 2.5;
    V3 toT = cam.target - ro;
    double tc = dot(toT, rd), d2 = dot(toT, toT) - tc * tc;
    if (d2 < BOUND * BOUND) {
        double half = std::sqrt(BOUND * BOUND - d2);
        double t = std::max(0.0, tc - half), tEnd = std::min(tc + half, tFloor);
        for (int i = 0; i < 160 && t < tEnd; ++i) {
            V3 p = ro + rd * t;
            double d = sdf(p, S);
            if (d < 0.0008 + 0.0007 * t) {
                -
                V3 n = normalAt(p, S);
                V3 vdir = rd * -1.0;
                double nv = clamp01(dot(n, vdir));
                double occ = 0.0, weight = 1.0;
                for (int j = 1; j <= 5; ++j) {                      
                    double h = 0.04 + 0.11 * j;
                    occ += (h - sdf(p + n * h, S)) * weight;
                    weight *= 0.8;
                }
                double ao = clamp01(1.0 - 1.6 * occ);
                double key = clamp01(dot(n, KEY)), fill = clamp01(dot(n, FILL));
                double spec = std::pow(clamp01(dot(n, unit(KEY + vdir))), 40.0);
                double L = 0.22 * ao + 0.52 * key * (0.4 + 0.6 * ao) + 0.20 * fill * ao + 0.20 * spec * key;
                L *= 0.55 + 0.45 * std::sqrt(nv);                    
                return {static_cast<float>(std::min(L, 0.93)), true};
            }
            t += d * 0.9;
        }
    }
    if (tFloor < 1e29) {                                              
        V3 f = ro + rd * tFloor;
        double sh = softShadow(f, KEY, S);
        double contact = std::exp(-4.0 * std::max(sdf(f, S), 0.0));
        double L = 1.0 - 0.32 * (1.0 - sh) - 0.22 * contact;
        return {static_cast<float>(std::max(L, 0.5)), false};
    }
    return {1.0f, false};
}

static std::vector<uint8_t> renderFrame(const Shape& S, const Cam& cam, int W, int H) {
    const double tanHalf = std::tan(0.5 * 30.0 * PI / 180.0);
    std::vector<float> L(static_cast<std::size_t>(W) * H);
    std::vector<uint8_t> hit(L.size());
    auto sample = [&](double px, double py) {
        return trace(cam, ((px / W) * 2.0 - 1.0) * tanHalf, (1.0 - (py / H) * 2.0) * tanHalf, S);
    };
    for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x) {
            Px s = sample(x + 0.5, y + 0.5);
            L[static_cast<std::size_t>(y) * W + x] = s.L;
            hit[static_cast<std::size_t>(y) * W + x] = s.hit ? 1 : 0;
        }
    
    for (int y = 1; y < H - 1; ++y)
        for (int x = 1; x < W - 1; ++x) {
            std::size_t i = static_cast<std::size_t>(y) * W + x;
            if (hit[i] == hit[i - 1] && hit[i] == hit[i + 1] && hit[i] == hit[i - W] && hit[i] == hit[i + W]) continue;
            double sum = L[i];
            for (int q = 0; q < 4; ++q) sum += sample(x + 0.5 + ((q & 1) ? 0.25 : -0.25), y + 0.5 + ((q & 2) ? 0.25 : -0.25)).L;
            L[i] = static_cast<float>(sum / 5.0);
        }
    std::vector<uint8_t> idx(L.size());
    for (std::size_t i = 0; i < L.size(); ++i) idx[i] = static_cast<uint8_t>(std::lround(clamp01(L[i]) * 255.0));
    return idx;
}


static void makePalette(double t, uint8_t pal[256][3]) {
    static const double theme[3][3][3] = {                    
        {{70, 36, 14}, {212, 146, 76}, {255, 233, 192}},
        {{74, 22, 56}, {206, 106, 138}, {255, 220, 226}},
        {{12, 34, 70}, {64, 136, 198}, {214, 236, 252}}};
    const double bg[3] = {248, 246, 242};
    double m = smooth(0.15, 0.85, t) * 2.0;                    
    int a = m < 1.0 ? 0 : 1;
    double f = m - a;
    for (int i = 0; i < 256; ++i) {
        double L = i / 255.0;
        for (int c = 0; c < 3; ++c) {
            double dark = mix(theme[a][0][c], theme[a + 1][0][c], f);
            double mid = mix(theme[a][1][c], theme[a + 1][1][c], f);
            double light = mix(theme[a][2][c], theme[a + 1][2][c], f);
            double v;
            if (L < 0.45) v = mix(dark, mid, L / 0.45);
            else if (L < 0.90) v = mix(mid, light, (L - 0.45) / 0.45);
            else v = mix(light, bg[c], (L - 0.90) / 0.10);
            pal[i][c] = static_cast<uint8_t>(std::lround(std::min(std::max(v, 0.0), 255.0)));
        }
    }
}


struct BitSink {
    std::vector<uint8_t> bytes;
    uint32_t acc = 0;
    int nbits = 0;
    void put(int code, int width) {
        acc |= static_cast<uint32_t>(code) << nbits;
        nbits += width;
        while (nbits >= 8) { bytes.push_back(static_cast<uint8_t>(acc & 0xFF)); acc >>= 8; nbits -= 8; }
    }
    void flush() { if (nbits > 0) { bytes.push_back(static_cast<uint8_t>(acc & 0xFF)); acc = 0; nbits = 0; } }
};

static std::vector<uint8_t> lzwEncode(const std::vector<uint8_t>& px) {
    const int CLEAR = 256, END = 257;
    std::vector<int16_t> table(4096 * 256, -1);                 
    BitSink out;
    int next = END + 1, width = 9;
    out.put(CLEAR, width);
    int prefix = px[0];
    for (std::size_t i = 1; i < px.size(); ++i) {
        int c = px[i];
        int16_t& slot = table[(static_cast<std::size_t>(prefix) << 8) | static_cast<std::size_t>(c)];
        if (slot >= 0) { prefix = slot; continue; }
        out.put(prefix, width);
        if (next < 4096) {
            slot = static_cast<int16_t>(next);
            if (next == (1 << width) && width < 12) ++width;
            ++next;
        } else {                                                
            out.put(CLEAR, width);
            std::fill(table.begin(), table.end(), static_cast<int16_t>(-1));
            next = END + 1;
            width = 9;
        }
        prefix = c;
    }
    out.put(prefix, width);
    if (next == (1 << width) && width < 12) ++width;
    out.put(END, width);
    out.flush();
    return out.bytes;
}

struct GifFrame { std::vector<uint8_t> data; uint8_t pal[256][3]; };

static void put16(std::ostream& o, int v) { o.put(static_cast<char>(v & 255)); o.put(static_cast<char>((v >> 8) & 255)); }

static bool writeGif(const char* file, int W, int H, const std::vector<GifFrame>& frames, const std::vector<std::pair<int, int>>& order) {
    std::ofstream o(file, std::ios::binary);
    if (!o) return false;
    o.write("GIF89a", 6);
    put16(o, W); put16(o, H);
    o.put(static_cast<char>(0x70)); o.put(0); o.put(0);         
    static const uint8_t loop[] = {0x21, 0xFF, 0x0B, 'N', 'E', 'T', 'S', 'C', 'A', 'P', 'E', '2', '.', '0', 0x03, 0x01, 0x00, 0x00, 0x00};
    o.write(reinterpret_cast<const char*>(loop), sizeof loop);  
    for (const auto& step : order) {
        const GifFrame& f = frames[static_cast<std::size_t>(step.first)];
        o.put(0x21); o.put(static_cast<char>(0xF9)); o.put(4); o.put(0x04);
        put16(o, step.second);                                  
        o.put(0); o.put(0);
        o.put(0x2C); put16(o, 0); put16(o, 0); put16(o, W); put16(o, H);
        o.put(static_cast<char>(0x87));                         
        o.write(reinterpret_cast<const char*>(f.pal), 768);
        o.put(8);                                               
        for (std::size_t i = 0; i < f.data.size(); i += 255) {
            std::size_t n = std::min<std::size_t>(255, f.data.size() - i);
            o.put(static_cast<char>(n));
            o.write(reinterpret_cast<const char*>(&f.data[i]), static_cast<std::streamsize>(n));
        }
        o.put(0);
    }
    o.put(0x3B);
    return static_cast<bool>(o);
}

int main() {
    const int W = 360, H = 360;         
    const int N = 41;                    

    std::vector<GifFrame> frames(N);
    std::ostringstream table;
    table << " stage     t    chi   b0  b1  b2   verdict\n";
    bool allSolidTori = true;

    std::cout << "Donut -> mug: building " << N << " stages and checking the topology of each...\n";
    for (int i = 0; i < N; ++i) {
        double t = static_cast<double>(i) / (N - 1);
        Shape S = shapeAt(t);

        Betti b = measureTopology(S);
        bool ok = (b.b0 == 1 && b.b1 == 1 && b.b2 == 0);
        allSolidTori = allSolidTori && ok;
        if (i % 4 == 0 || !ok) {
            table << std::setw(5) << i << std::setw(8) << std::fixed << std::setprecision(2) << t << std::setw(7) << b.chi
                  << std::setw(6) << b.b0 << std::setw(4) << b.b1 << std::setw(4) << b.b2
                  << (ok ? "   one solid with one hole (donut-type)" : "   TOPOLOGY CHANGED") << "\n";
        }

        frames[static_cast<std::size_t>(i)].data = lzwEncode(renderFrame(S, cameraAt(t), W, H));
        makePalette(t, frames[static_cast<std::size_t>(i)].pal);
        std::cout << "\r  stage " << (i + 1) << "/" << N << std::flush;
    }
    std::cout << "\n\n" << table.str() << "\n";
    std::cout << (allSolidTori ? "All stages have (b0,b1,b2) = (1,1,0): the hole is never closed or added.\n"
                               : "WARNING: the topology changed at some stage.\n");

    std::vector<std::pair<int, int>> order;
    for (int i = 0; i < N; ++i) order.push_back({i, (i == 0 || i == N - 1) ? 90 : 7});
    for (int i = N - 2; i >= 1; --i) order.push_back({i, 7});

    if (!writeGif("donut_to_mug.gif", W, H, frames, order)) {
        std::cerr << "Cannot write donut_to_mug.gif\n";
        return 1;
    }
    std::ofstream html("donut_to_mug.html");
    html << "<!doctype html>\n<html><head><meta charset=\"utf-8\"><title>Donut to mug</title>\n"
         << "<style>body{font-family:'Segoe UI',Arial,sans-serif;text-align:center;background:#f4f6f8;color:#222}"
         << "img{border:1px solid #ccd;border-radius:8px;background:#fff}pre{display:inline-block;text-align:left;background:#fff;"
         << "border:1px solid #ccd;border-radius:8px;padding:10px 16px}</style></head><body>\n"
         << "<h2>Donut &rarr; mug: no tearing, no gluing, no new holes</h2>\n"
         << "<img src=\"donut_to_mug.gif\" width=\"" << W << "\" height=\"" << H << "\"><br>\n"
         << "<p>Topology check of the stages (b0 = pieces, b1 = tunnels, b2 = hollow cavities):</p>\n"
         << "<pre>" << table.str() << "</pre>\n</body></html>\n";
    html.close();

    std::error_code ec;
    const auto full = std::filesystem::absolute("donut_to_mug.html", ec);
    std::cout << "Wrote " << (ec ? std::string("donut_to_mug.html") : full.string()) << "\n";
#ifdef _WIN32
    std::system("start \"\" \"donut_to_mug.html\"");           
#else
    std::cout << "Open it in a web browser to watch the animation.\n";
#endif
    return allSolidTori ? 0 : 1;
}
