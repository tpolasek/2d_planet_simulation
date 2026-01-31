#include "simulation.h"
#include "utility.h"
#include <cmath>
#include <random>
#include <algorithm>
#include <iostream>

int Simulation::s_quadtreeNodeCount = 0;

Simulation::Simulation(int width, int height, int num_particles, double elasticity)
    : width(width), height(height), num_particles(num_particles), elasticity(elasticity) {
    occupancy.resize(width * height, -1);
    
    // Initialize uniform grid
    cellSize = 32;
    gridCols = (width + cellSize - 1) / cellSize;
    gridRows = (height + cellSize - 1) / cellSize;
    grid.resize(gridCols * gridRows);
    
    spawnParticles();
}

void Simulation::spawnParticles() {
    std::random_device rd;
    std::mt19937 gen(rd());
    
    // Spawn within a 200x200 box at the center
    std::uniform_real_distribution<> dis_x(100, width - 100);
    std::uniform_real_distribution<> dis_y(100, height - 100);
    std::uniform_real_distribution<> dis_v(-0.5, 0.5);

    // occupancy grid is already sized width*height and filled with -1 (by constructor)
    // We'll use it to ensure each pixel holds at most one particle
    for (int i = 0; i < num_particles; ++i) {
        int attempts = 0;
        const int max_attempts = 100;
        double x, y;
        int px, py, idx;
        bool placed = false;
        
        while (!placed && attempts < max_attempts) {
            x = dis_x(gen);
            y = dis_y(gen);
            px = std::clamp((int)std::round(x), 0, width - 1);
            py = std::clamp((int)std::round(y), 0, height - 1);
            idx = py * width + px;
            
            if (occupancy[idx] == -1) {
                placed = true;
                occupancy[idx] = i; // mark pixel occupied by this particle index
            } else {
                ++attempts;
            }
        }
        
        // If we failed to find an empty pixel after many attempts, just place anyway
        // (overlap will be resolved later by collision detection)
        double centerX = width * 0.5;
        double centerY = height * 0.5;
        double dx = x - centerX;
        double dy = y - centerY;
        double dist_sq = dx*dx + dy*dy;
        double vx, vy;
        if (dist_sq > 1e-12) {
            // clockwise tangent: (dy, -dx)
            double tx = dy;
            double ty = -dx;
            double inv_len = 1.0 / std::sqrt(dist_sq);
            double speed = 0.3; // constant tangential speed
            vx = tx * inv_len * speed;
            vy = ty * inv_len * speed;
        } else {
            // particle at center: give random direction
            vx = dis_v(gen);
            vy = dis_v(gen);
        }

        particles.push_back({
            x, y,
            x - vx, y - vy, // prev_x, prev_y for initial velocity
            0.0, 0.0,
            1.0,
            255, 0, 0
        });
    }
}

double FastInvSqrt(double x) {
    double xhalf = 0.5 * x;
    long long i = *(long long*)&x;         // Treat bits as integer
    i = 0x5fe6eb50c7b537a9LL - (i >> 1);   // The magic constant for double
    x = *(double*)&i;                      // Treat bits back as double
    x = x * (1.5 - xhalf * x * x);         // Newton-Raphson iteration
    // x = x * (1.5 - xhalf * x * x);      // Optional: second iteration
    return x;
}

void Simulation::applyGravity() {
    buildQuadtree();
    applyGravityUsingQuadtree();
}

void Simulation::verletStep(double dt) {
    const size_t n = particles.size();
    const double dt2 = dt * dt;

    //#pragma omp parallel
    for (size_t i = 0; i < n; ++i) {
        Particle& p = particles[i];
        double temp_x = p.x;
        double temp_y = p.y;

        double new_x = 2.0 * p.x - p.prev_x + p.ax * dt2;
        double new_y = 2.0 * p.y - p.prev_y + p.ay * dt2;

        // Validate finite positions; if invalid, revert to previous position and zero velocity
        if (!std::isfinite(new_x) || !std::isfinite(new_y)) {
            new_x = p.x;
            new_y = p.y;
            p.prev_x = p.x;
            p.prev_y = p.y;
        }

        p.x = new_x;
        p.y = new_y;

        p.prev_x = temp_x;
        p.prev_y = temp_y;

        p.ax = 0.0;
        p.ay = 0.0;

        // Boundary checks
        if (p.x < 0) { p.x = 0; p.prev_x = p.x + (p.x - p.prev_x) * -0.5; }
        if (p.x >= width) { p.x = width - 1; p.prev_x = p.x + (p.x - p.prev_x) * -0.5; }
        if (p.y < 0) { p.y = 0; p.prev_y = p.y + (p.y - p.prev_y) * -0.5; }
        if (p.y >= height) { p.y = height - 1; p.prev_y = p.y + (p.y - p.prev_y) * -0.5; }
    }
}

void Simulation::resolvePixelCollisions() {
    // Clear occupancy grid
    std::fill(occupancy.begin(), occupancy.end(), -1);

    const double e = elasticity;
    const size_t n = particles.size();

    for (size_t i = 0; i < n; ++i) {
        int px = std::clamp((int)std::round(particles[i].x), 0, width - 1);
        int py = std::clamp((int)std::round(particles[i].y), 0, height - 1);
        int idx = py * width + px;

        if (occupancy[idx] != -1) {
            // Collision with particle j
            size_t j = occupancy[idx];
            
            // Load particle data into locals
            double xi = particles[i].x, yi = particles[i].y;
            double xj = particles[j].x, yj = particles[j].y;
            double dx = xj - xi;
            double dy = yj - yi;
            double dist_sq = dx * dx + dy * dy + 1e-12; // avoid zero
            if (!std::isfinite(dist_sq) || dist_sq < 1e-24) {
                // Particles are effectively at same position; skip collision
                occupancy[idx] = i;
                continue;
            }
            double inv_dist = 1.0 / std::sqrt(dist_sq);
            double nx = dx * inv_dist;
            double ny = dy * inv_dist;

            // Velocities (Verlet)
            double v1x = xi - particles[i].prev_x;
            double v1y = yi - particles[i].prev_y;
            double v2x = xj - particles[j].prev_x;
            double v2y = yj - particles[j].prev_y;

            // Normal velocities
            double v1n = v1x * nx + v1y * ny;
            double v2n = v2x * nx + v2y * ny;

            // Masses
            double m1 = particles[i].mass;
            double m2 = particles[j].mass;
            double total_mass = m1 + m2;
            double inv_total_mass = 1.0 / total_mass;

            // Elastic collision response
            double v1n_after = ((m1 - e * m2) * v1n + (1.0 + e) * m2 * v2n) * inv_total_mass;
            double v2n_after = ((m2 - e * m1) * v2n + (1.0 + e) * m1 * v1n) * inv_total_mass;

            // Velocity deltas
            double dv1n = v1n_after - v1n;
            double dv2n = v2n_after - v2n;

            // Update previous positions (which encode velocity)
            particles[i].prev_x = xi - (v1x + dv1n * nx);
            particles[i].prev_y = yi - (v1y + dv1n * ny);
            particles[j].prev_x = xj - (v2x + dv2n * nx);
            particles[j].prev_y = yj - (v2y + dv2n * ny);

            // Separate positions by half pixel along normal
            const double separation = 0.5;
            double sx = separation * nx;
            double sy = separation * ny;
            particles[i].x = xi - sx;
            particles[i].y = yi - sy;
            particles[i].prev_x -= sx;
            particles[i].prev_y -= sy;
            particles[j].x = xj + sx;
            particles[j].y = yj + sy;
            particles[j].prev_x += sx;
            particles[j].prev_y += sy;
        } else {
            occupancy[idx] = i;
        }
    }
}

void Simulation::render(SDL_Renderer* renderer) {
    const double maxSpeed = 1.0; // expected maximum speed for hue scaling
    const double saturation = 1.0;
    const double value = 1.0;
    for (const auto& p : particles) {
        double vx = p.x - p.prev_x;
        double vy = p.y - p.prev_y;
        double speed = std::sqrt(vx*vx + vy*vy);
        double hue = (speed / maxSpeed) * 360.0;
        // clamp hue to [0,360)
        if (hue < 0.0) hue = 0.0;
        if (hue >= 360.0) hue = 359.999;
        uint8_t r, g, b;
        HSVtoRGB(hue, saturation, value, r, g, b);
        SDL_SetRenderDrawColor(renderer, r, g, b, 255);
        SDL_RenderDrawPoint(renderer, (int)std::round(p.x), (int)std::round(p.y));
    }
}

void Simulation::clearUniformGrid() {
    for (auto& cell : grid) {
        cell.clear();
    }
}

void Simulation::updateUniformGrid() {
    clearUniformGrid();
    for (size_t i = 0; i < particles.size(); ++i) {
        int cellX = static_cast<int>(particles[i].x) / cellSize;
        int cellY = static_cast<int>(particles[i].y) / cellSize;
        if (cellX >= 0 && cellX < gridCols && cellY >= 0 && cellY < gridRows) {
            grid[cellY * gridCols + cellX].push_back(i);
        }
    }
}

// QuadtreeNode constructor
Simulation::QuadtreeNode::QuadtreeNode(double cx, double cy, double sz)
    : centerX(cx), centerY(cy), size(sz), totalMass(0.0), centerOfMassX(0.0), centerOfMassY(0.0),
      particleIndex(-1), particleX(0.0), particleY(0.0), particleMass(0.0) {
    Simulation::s_quadtreeNodeCount++;
    for (int i = 0; i < 4; ++i) children[i].reset();
}

Simulation::QuadtreeNode::~QuadtreeNode() {
    Simulation::s_quadtreeNodeCount--;
}

void Simulation::QuadtreeNode::insert(int pIdx, double px, double py, double pmass, int depth) {
    // Safety check: if depth exceeds maximum, treat as leaf and overwrite existing particle
    if (depth >= MAX_DEPTH) {
        particleIndex = pIdx;
        particleX = px;
        particleY = py;
        particleMass = pmass;
        return;
    }
    
    // If this node is empty and has no children, store particle here
    if (particleIndex == -1 && !children[0]) {
        particleIndex = pIdx;
        particleX = px;
        particleY = py;
        particleMass = pmass;
        return;
    }
    
    // If this node has children, insert into appropriate child
    if (children[0]) {
        int quadrant = 0;
        if (px >= centerX) quadrant |= 1; // east
        if (py >= centerY) quadrant |= 2; // south
        children[quadrant]->insert(pIdx, px, py, pmass, depth + 1);
        return;
    }
    
    // Otherwise, node contains a particle, need to subdivide
    // Store existing particle data
    int existingIdx = particleIndex;
    double existingX = particleX;
    double existingY = particleY;
    double existingMass = particleMass;
    particleIndex = -1;
    
    // Create four children
    double halfSize = size * 0.5;
    double childSize = halfSize;
    double left = centerX - halfSize;
    double right = centerX + halfSize;
    double top = centerY - halfSize;
    double bottom = centerY + halfSize;
    
    children[0].reset(new QuadtreeNode(left + childSize*0.5, top + childSize*0.5, childSize));   // NW
    children[1].reset(new QuadtreeNode(right - childSize*0.5, top + childSize*0.5, childSize));  // NE
    children[2].reset(new QuadtreeNode(left + childSize*0.5, bottom - childSize*0.5, childSize)); // SW
    children[3].reset(new QuadtreeNode(right - childSize*0.5, bottom - childSize*0.5, childSize)); // SE
    
    // Re-insert existing particle
    int quadrant = 0;
    if (existingX >= centerX) quadrant |= 1;
    if (existingY >= centerY) quadrant |= 2;
    children[quadrant]->insert(existingIdx, existingX, existingY, existingMass, depth + 1);
    
    // Insert new particle
    quadrant = 0;
    if (px >= centerX) quadrant |= 1;
    if (py >= centerY) quadrant |= 2;
    children[quadrant]->insert(pIdx, px, py, pmass, depth + 1);
}

void Simulation::QuadtreeNode::computeMassDistribution() {
    if (particleIndex != -1 && !children[0]) {
        // Leaf with a single particle
        totalMass = particleMass;
        centerOfMassX = particleX;
        centerOfMassY = particleY;
        return;
    }
    
    // Internal node: compute from children
    totalMass = 0.0;
    double weightedX = 0.0;
    double weightedY = 0.0;
    for (int i = 0; i < 4; ++i) {
        if (children[i]) {
            children[i]->computeMassDistribution();
            totalMass += children[i]->totalMass;
            weightedX += children[i]->totalMass * children[i]->centerOfMassX;
            weightedY += children[i]->totalMass * children[i]->centerOfMassY;
        }
    }
    if (totalMass > 0.0) {
        centerOfMassX = weightedX / totalMass;
        centerOfMassY = weightedY / totalMass;
    }
}

void Simulation::QuadtreeNode::applyGravityToParticle(int pIdx, double px, double py, double pmass,
                                                      double& ax, double& ay, double theta, double G, double EPS) const {
    // Skip self-interaction if this leaf contains the same particle
    if (particleIndex == pIdx && !children[0]) return;
    
    double dx = centerOfMassX - px;
    double dy = centerOfMassY - py;
    double dist_sq = dx * dx + dy * dy + EPS * EPS;
    double dist = std::sqrt(dist_sq);
    
    // Opening criterion
    if (size / dist < theta || (!children[0] && particleIndex != -1)) {
        // Treat as a single mass
        double inv_dist_sq = 1.0 / dist_sq;
        double force_mag = G * pmass * totalMass * inv_dist_sq;
        double factor = force_mag / dist;
        ax += factor * dx;
        ay += factor * dy;
    } else {
        // Recursively process children
        for (int i = 0; i < 4; ++i) {
            if (children[i]) {
                children[i]->applyGravityToParticle(pIdx, px, py, pmass, ax, ay, theta, G, EPS);
            }
        }
    }
}

void Simulation::buildQuadtree() {
    // Determine root bounds covering all particles, ignoring non‑finite positions
    double minX = width, maxX = 0, minY = height, maxY = 0;
    bool hasValid = false;
    for (const auto& p : particles) {
        if (!std::isfinite(p.x) || !std::isfinite(p.y)) continue;
        hasValid = true;
        if (p.x < minX) minX = p.x;
        if (p.x > maxX) maxX = p.x;
        if (p.y < minY) minY = p.y;
        if (p.y > maxY) maxY = p.y;
    }
    // If no valid particles, create a dummy root covering the whole screen
    if (!hasValid) {
        minX = 0; maxX = width;
        minY = 0; maxY = height;
    }
    // Add a little padding
    double padding = 1.0;
    minX -= padding; maxX += padding;
    minY -= padding; maxY += padding;
    double size = std::max(maxX - minX, maxY - minY);
    // Ensure size is finite and positive
    if (!std::isfinite(size) || size <= 0.0) {
        size = width;
    }
    double centerX = (minX + maxX) * 0.5;
    double centerY = (minY + maxY) * 0.5;
    if (!std::isfinite(centerX) || !std::isfinite(centerY)) {
        centerX = width * 0.5;
        centerY = height * 0.5;
    }
    
    quadtreeRoot.reset(new QuadtreeNode(centerX, centerY, size));
    
    // Insert all particles with finite positions
    for (size_t i = 0; i < particles.size(); ++i) {
        const Particle& p = particles[i];
        if (std::isfinite(p.x) && std::isfinite(p.y)) {
            quadtreeRoot->insert(i, p.x, p.y, p.mass);
        }
    }
    
    // Compute mass distribution
    quadtreeRoot->computeMassDistribution();
}

void Simulation::applyGravityUsingQuadtree() {
    if (!quadtreeRoot) return;
    
    const double G_local = G;
    const double EPS_local = EPS;
    const double theta_local = theta;
    
    for (size_t i = 0; i < particles.size(); ++i) {
        Particle& p = particles[i];
        double ax = 0.0, ay = 0.0;
        quadtreeRoot->applyGravityToParticle(i, p.x, p.y, p.mass, ax, ay, theta_local, G_local, EPS_local);
        p.ax += ax;
        p.ay += ay;
    }
}

void Simulation::logMemoryStats() const {
    size_t occupancy_capacity = occupancy.capacity() * sizeof(int);
    size_t particles_capacity = particles.capacity() * sizeof(Particle);
    size_t grid_capacity = grid.capacity() * sizeof(std::vector<int>);
    size_t grid_inner_capacity = 0;
    for (const auto& cell : grid) {
        grid_inner_capacity += cell.capacity() * sizeof(int);
    }
    std::cerr << "Memory stats: "
              << "occupancy capacity=" << occupancy_capacity << "B, "
              << "particles capacity=" << particles_capacity << "B, "
              << "grid outer capacity=" << grid_capacity << "B, "
              << "grid inner capacity=" << grid_inner_capacity << "B"
              << std::endl;
}

double Simulation::getMaxSpeed() const {
    double maxSpeedSq = 0.0;
    for (const auto& p : particles) {
        double vx = p.x - p.prev_x;
        double vy = p.y - p.prev_y;
        double speedSq = vx*vx + vy*vy;
        if (speedSq > maxSpeedSq) maxSpeedSq = speedSq;
    }
    return std::sqrt(maxSpeedSq);
}
