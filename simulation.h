#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include <SDL2/SDL.h>
#include "particle.h"
#include <memory>
#include <cmath>

class Simulation {
public:
    Simulation(int width, int height, int num_particles, double elasticity = 1.0);
    void spawnParticles();
    void applyGravity();
    void verletStep(double dt);
    void resolvePixelCollisions();
    void render(SDL_Renderer* renderer);
    static int getQuadtreeNodeCount() { return s_quadtreeNodeCount; }
    void logMemoryStats() const;
    double getMaxSpeed() const;

private:
    int width, height;
    int num_particles;
    std::vector<Particle> particles;
    std::vector<int> occupancy; // width * height grid
    double elasticity;

    const double G = 100.0; // Gravitational constant (tuned for visual effect)
    const double EPS = 1.0; // Softening factor

    // Uniform grid for broad-phase collision detection
    int cellSize = 32;
    int gridCols, gridRows;
    std::vector<std::vector<int>> grid; // grid[cellY * gridCols + cellX] = list of particle indices
    void updateUniformGrid();
    void clearUniformGrid();

    // Barnes-Hut quadtree for gravity
    struct QuadtreeNode {
        static constexpr int MAX_DEPTH = 10000;
        double centerX, centerY; // center of this node
        double size; // width and height of this node (square)
        double totalMass;
        double centerOfMassX, centerOfMassY;
        int particleIndex = -1; // if leaf with a single particle
        double particleX, particleY, particleMass; // stored copy for subdivision
        std::unique_ptr<QuadtreeNode> children[4]; // NW, NE, SW, SE

        QuadtreeNode(double cx, double cy, double sz);
        ~QuadtreeNode();
        void insert(int pIdx, double px, double py, double pmass, int depth = 0);
        void computeMassDistribution();
        void applyGravityToParticle(int pIdx, double px, double py, double pmass,
                                    double& ax, double& ay, double theta, double G, double EPS) const;
    };
    std::unique_ptr<QuadtreeNode> quadtreeRoot;
    const double theta = 0.5; // Barnes-Hut opening angle
    void buildQuadtree();
    void applyGravityUsingQuadtree();

    static int s_quadtreeNodeCount;
};

#endif // SIMULATION_H
