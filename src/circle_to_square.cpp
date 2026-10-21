#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

struct Pt { double x, y; };

static const double PI = 3.14159265358979323846;

Pt circleToSquare(Pt p) {
    double m = std::max(std::fabs(p.x), std::fabs(p.y));   
    if (m == 0.0) return {0.0, 0.0};
    double k = std::hypot(p.x, p.y) / m;                   
    return {p.x * k, p.y * k};
}

Pt squareToCircle(Pt q) {
    double r = std::hypot(q.x, q.y);
    if (r == 0.0) return {0.0, 0.0};
    double k = std::max(std::fabs(q.x), std::fabs(q.y)) / r;
    return {q.x * k, q.y * k};
}


Pt morph(Pt p, double t) {
    Pt q = circleToSquare(p);
    return {(1.0 - t) * p.x + t * q.x, (1.0 - t) * p.y + t * q.y};
}


Pt onCircle(double angle, double radius = 1.0) {
    return {radius * std::cos(angle), radius * std::sin(angle)};
}


std::vector<Pt> ring(double radius, int n) {
    std::vector<Pt> v;
    for (int k = 0; k < n; ++k) v.push_back(onCircle(2.0 * PI * k / n, radius));
    return v;
}

Pt onSquare(double s) {
    int side = static_cast<int>(s);
    double u = 2.0 * (s - side) - 1.0;                    
    switch (side) {
        case 0:  return { 1.0,   u};
        case 1:  return {  -u, 1.0};
        case 2:  return {-1.0,  -u};
        default: return {   u, -1.0};
    }
}


std::vector<Pt> mapAll(const std::vector<Pt>& v) {
    std::vector<Pt> w;
    for (const Pt& p : v) w.push_back(circleToSquare(p));
    return w;
}


double dist(Pt a, Pt b) { return std::hypot(a.x - b.x, a.y - b.y); }


bool windsOnce(const std::vector<Pt>& poly) {
    double total = 0.0;
    for (std::size_t i = 0; i < poly.size(); ++i) {
        Pt a = poly[i], b = poly[(i + 1) % poly.size()];
        double turn = std::atan2(a.x * b.y - a.y * b.x, a.x * b.x + a.y * b.y); 
        if (turn <= 0.0) return false;
        total += turn;
    }
    return std::fabs(total - 2.0 * PI) < 1e-9;
}

std::string pointsAttr(const std::vector<Pt>& v) {
    std::ostringstream s;
    s << std::fixed << std::setprecision(4);
    for (const Pt& p : v) s << p.x << ',' << -p.y << ' '; 
    return s.str();
}


void animatedShape(std::ostream& out, const char* tag, const char* style,
                   const std::vector<Pt>& A, const std::vector<Pt>& B) {
    const std::string a = pointsAttr(A), b = pointsAttr(B);
    out << "<" << tag << " " << style << " points=\"" << a << "\">\n"
        << "  <animate attributeName=\"points\" dur=\"7s\" repeatCount=\"indefinite\"\n"
        << "    calcMode=\"spline\" keyTimes=\"0;0.15;0.5;0.65;1\"\n"
        << "    keySplines=\"0 0 1 1;0.6 0 0.4 1;0 0 1 1;0.6 0 0.4 1\"\n"
        << "    values=\"" << a << ";" << a << ";" << b << ";" << b << ";" << a << "\"/>\n"
        << "</" << tag << ">\n";
}

int main() {
   
    const int M = 100000;
    double offSquare = 0.0, errCircle = 0.0, errSquare = 0.0;
    for (int i = 0; i < M; ++i) {
        Pt p = onCircle(2.0 * PI * i / M);
        Pt q = circleToSquare(p);
        offSquare = std::max(offSquare, std::fabs(std::max(std::fabs(q.x), std::fabs(q.y)) - 1.0));
        errCircle = std::max(errCircle, dist(squareToCircle(q), p));

        Pt s = onSquare(4.0 * i / M);
        errSquare = std::max(errSquare, dist(circleToSquare(squareToCircle(s)), s));
    }

    const int N = 360;                                    
    const std::vector<Pt> circle = ring(1.0, N);
    bool simpleAlways = true;                             
    for (int i = 0; i <= 10; ++i) {
        std::vector<Pt> shape;
        for (const Pt& p : circle) shape.push_back(morph(p, i / 10.0));
        simpleAlways = simpleAlways && windsOnce(shape);
    }

    const bool ok = offSquare < 1e-12 && errCircle < 1e-12 && errSquare < 1e-12 && simpleAlways;

    std::cout << "Circle -> square with a homeomorphism\n"
              << "  f(x,y) = (x,y) * sqrt(x^2+y^2) / max(|x|,|y|)    circle -> square\n"
              << "  g(u,v) = (u,v) * max(|u|,|v|) / sqrt(u^2+v^2)    square -> circle\n\n"
              << std::scientific << std::setprecision(1)
              << "Checks (" << M << " sample points each):\n"
              << "  f(circle) lies on the square's outline : max error " << offSquare << "\n"
              << "  g(f(p)) = p  for p on the circle       : max error " << errCircle << "\n"
              << "  f(g(q)) = q  for q on the square       : max error " << errSquare << "\n"
              << "  one simple closed curve at t = 0, 0.1, ..., 1 : " << (simpleAlways ? "yes" : "NO") << "\n"
              << (ok ? "=> all checks passed\n\n" : "=> SOMETHING IS WRONG\n\n");

    const char* file = "circle_to_square.html";
    std::ofstream out(file);
    if (!out) {
        std::cerr << "Cannot write " << file << "\n";
        return 1;
    }

    out << "<!doctype html>\n<html><head><meta charset=\"utf-8\"><title>Circle to square</title>\n"
        << "<style>body{font-family:'Segoe UI',Arial,sans-serif;text-align:center;background:#f4f6f8;color:#222}"
        << "svg{background:#fff;border:1px solid #ccd;border-radius:8px;width:min(90vw,80vh);height:auto}</style>\n"
        << "</head><body>\n<h2>Circle &rarr; square: no tearing, no gluing</h2>\n"
        << "<svg viewBox=\"-1.25 -1.25 2.5 2.5\" xmlns=\"http://www.w3.org/2000/svg\">\n";

    for (double t : {0.2, 0.4, 0.6, 0.8}) {
        std::vector<Pt> shape;
        for (const Pt& p : circle) shape.push_back(morph(p, t));
        out << "<polygon fill=\"none\" stroke=\"#e08a1e\" stroke-opacity=\"0.55\" stroke-width=\"0.008\" points=\""
            << pointsAttr(shape) << "\"/>\n";
    }

    // The inside of the disk, carried along by the same map
    const char* thin = "fill=\"none\" stroke=\"#9db4d0\" stroke-width=\"0.006\"";
    for (double r : {0.25, 0.5, 0.75}) {
        std::vector<Pt> a = ring(r, N);
        animatedShape(out, "polygon", thin, a, mapAll(a));
    }
    
    std::vector<Pt> spokes;
    for (int k = 0; k < 24; ++k) {
        spokes.push_back({0.0, 0.0});
        spokes.push_back(onCircle(2.0 * PI * k / 24));
    }
    animatedShape(out, "polyline", thin, spokes, mapAll(spokes));

    // The outline itself, on top.
    animatedShape(out, "polygon",
                  "fill=\"#4682dc\" fill-opacity=\"0.12\" stroke=\"#1f5fbf\" stroke-width=\"0.02\" stroke-linejoin=\"round\"",
                  circle, mapAll(circle));

    out << "</svg>\n"
        << "<p>Blue: the shape while it deforms &middot; Orange: snapshots at t = 0.2, 0.4, 0.6, 0.8 &middot; "
        << "Grey: rings and spokes inside the disk, carried along by the same map.</p>\n"
        << "</body></html>\n";
    out.close();

    std::error_code ec;
    const auto full = std::filesystem::absolute(file, ec);
    std::cout << "Wrote " << (ec ? std::string(file) : full.string()) << "\n";
#ifdef _WIN32
    std::system("start \"\" \"circle_to_square.html\"");   
#else
    std::cout << "Open it in a web browser to watch the animation.\n";
#endif
    return ok ? 0 : 1;
}
