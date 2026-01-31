#include "utility.h"
#include <cmath>

void HSVtoRGB(double h, double s, double v, uint8_t& r, uint8_t& g, uint8_t& b) {
    // wrap hue to [0,360)
    while (h < 0.0) h += 360.0;
    while (h >= 360.0) h -= 360.0;
    double c = v * s;
    double x = c * (1.0 - std::fabs(std::fmod(h / 60.0, 2.0) - 1.0));
    double m = v - c;
    double r1, g1, b1;
    if (h < 60.0) { r1 = c; g1 = x; b1 = 0; }
    else if (h < 120.0) { r1 = x; g1 = c; b1 = 0; }
    else if (h < 180.0) { r1 = 0; g1 = c; b1 = x; }
    else if (h < 240.0) { r1 = 0; g1 = x; b1 = c; }
    else if (h < 300.0) { r1 = x; g1 = 0; b1 = c; }
    else { r1 = c; g1 = 0; b1 = x; }
    r = static_cast<uint8_t>((r1 + m) * 255.0);
    g = static_cast<uint8_t>((g1 + m) * 255.0);
    b = static_cast<uint8_t>((b1 + m) * 255.0);
}
