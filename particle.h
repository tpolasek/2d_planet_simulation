#ifndef PARTICLE_H
#define PARTICLE_H

#include <cstdint>

struct Particle {
    double x, y;          // current position
    double prev_x, prev_y;// previous position (for Verlet)
    double ax, ay;        // accumulated acceleration
    double mass;          // mass
    uint8_t r, g, b;      // color
};

#endif // PARTICLE_H
