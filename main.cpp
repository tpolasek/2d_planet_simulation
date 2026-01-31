#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include "simulation.h"

int main(int argc, char* argv[]) {
    int num_particles = 8000; // default

    if (argc > 1) {
        char* endptr;
        long n = std::strtol(argv[1], &endptr, 10);
        if (endptr == argv[1] || n <= 0) {
            std::cerr << "Usage: " << argv[0] << " <N>\n";
            std::cerr << "  N: number of particles (1-500)\n";
            return 1;
        }
        num_particles = static_cast<int>(n);
    }

    const int WIDTH = 1000;
    const int HEIGHT = 1000;
    const double DT = 0.016; // ~60 FPS

    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        std::cerr << "SDL_Init failed: " << SDL_GetError() << "\n";
        return 1;
    }

    if (TTF_Init() != 0) {
        std::cerr << "TTF_Init failed: " << TTF_GetError() << "\n";
        SDL_Quit();
        return 1;
    }

    SDL_Window* window = SDL_CreateWindow(
        "Pixel Gravity Simulation",
        SDL_WINDOWPOS_CENTERED,
        SDL_WINDOWPOS_CENTERED,
        WIDTH,
        HEIGHT,
        SDL_WINDOW_SHOWN
    );

    if (!window) {
        std::cerr << "SDL_CreateWindow failed: " << SDL_GetError() << "\n";
        SDL_Quit();
        return 1;
    }

    SDL_Renderer* renderer = SDL_CreateRenderer(
        window,
        -1,
        SDL_RENDERER_ACCELERATED
    );

    if (!renderer) {
        std::cerr << "SDL_CreateRenderer failed: " << SDL_GetError() << "\n";
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    Simulation sim(WIDTH, HEIGHT, num_particles, 0.5);

    // Load font for FPS display
    TTF_Font* font = TTF_OpenFont("/System/Library/Fonts/Helvetica.ttc", 16);
    if (!font) {
        std::cerr << "Failed to load font: " << TTF_GetError() << "\n";
        TTF_Quit();
        SDL_Quit();
        return 1;
    }

    bool running = true;
    Uint32 frame_start;
    Uint32 frame_time;
    Uint32 last_fps_update = 0;
    int frame_count = 0;
    double fps = 0.0;

    int frame_counter = 0;
    while (running) {
        frame_start = SDL_GetTicks();

        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                running = false;
            } else if (event.type == SDL_KEYDOWN) {
                if (event.key.keysym.sym == SDLK_ESCAPE) {
                    running = false;
                }
            }
        }

        sim.applyGravity();
        sim.verletStep(DT);
        sim.resolvePixelCollisions();

        // Log quadtree node count every 60 frames
        frame_counter++;
        if (frame_counter % 60 == 0) {
            int node_count = Simulation::getQuadtreeNodeCount();
            double maxSpeed = sim.getMaxSpeed();
            std::cerr << "Frame " << frame_counter << ": Quadtree nodes = " << node_count
                      << ", max speed = " << maxSpeed << std::endl;
            sim.logMemoryStats();
        }

        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        sim.render(renderer);

        // Calculate FPS
        Uint32 current_time = SDL_GetTicks();
        frame_count++;
        if (current_time - last_fps_update >= 500) { // Update every 500ms
            fps = frame_count * 1000.0 / (current_time - last_fps_update);
            frame_count = 0;
            last_fps_update = current_time;
        }

        // Render FPS text in top right corner
        std::ostringstream fps_text;
        fps_text << "FPS: " << static_cast<int>(fps);
        SDL_Color text_color = {255, 255, 255, 255};
        SDL_Surface* surface = TTF_RenderText_Solid(font, fps_text.str().c_str(), text_color);
        if (surface) {
            SDL_Texture* texture = SDL_CreateTextureFromSurface(renderer, surface);
            if (texture) {
                int text_width = surface->w;
                int text_height = surface->h;
                SDL_Rect dest_rect = {WIDTH - text_width - 10, 10, text_width, text_height};
                SDL_RenderCopy(renderer, texture, NULL, &dest_rect);
                SDL_DestroyTexture(texture);
            }
            SDL_FreeSurface(surface);
        }

        SDL_RenderPresent(renderer);

        frame_time = SDL_GetTicks() - frame_start;
        if (frame_time < 16) {
            SDL_Delay(16 - frame_time);
        }
    }

    TTF_CloseFont(font);
    TTF_Quit();
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
