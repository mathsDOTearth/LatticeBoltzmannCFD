#include <SDL2/SDL.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "lbm.h"

#define WINDOW_WIDTH 800
#define WINDOW_HEIGHT 800

int main() {
    // Initialise LBM
    LBM lbm;
    int nx = 200, ny = 200;  // Increased grid resolution
    initialise(&lbm, nx, ny);

    // Initialise SDL
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        fprintf(stderr, "SDL could not initialize! SDL_Error: %s\n", SDL_GetError());
        return 1;
    }

    SDL_Window *window = SDL_CreateWindow("LBM Simulation Visualisation",
                                          SDL_WINDOWPOS_UNDEFINED,
                                          SDL_WINDOWPOS_UNDEFINED,
                                          WINDOW_WIDTH, WINDOW_HEIGHT,
                                          SDL_WINDOW_SHOWN);

    if (window == NULL) {
        fprintf(stderr, "Window could not be created! SDL_Error: %s\n", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    SDL_Renderer *renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

    if (renderer == NULL) {
        fprintf(stderr, "Renderer could not be created! SDL_Error: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    SDL_Event event;
    int running = 1;
    int time_steps = 1000;

    for (int t = 0; t < time_steps && running; t++) {
        // Process events
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                running = 0;
            }
        }

        // Update the simulation
        compute_macroscopic_variables(&lbm);
        collision_step(&lbm);
        streaming_step(&lbm);
        boundary_conditions(&lbm);

        // Compute min and max density values
        double min_rho = lbm.rho[0][0];
        double max_rho = lbm.rho[0][0];

        #pragma omp parallel for collapse(2) reduction(min:min_rho) reduction(max:max_rho)
        for (int x = 0; x < nx; x++) {
            for (int y = 0; y < ny; y++) {
                if (lbm.rho[x][y] < min_rho) min_rho = lbm.rho[x][y];
                if (lbm.rho[x][y] > max_rho) max_rho = lbm.rho[x][y];
            }
        }

        // Apply optional smoothing to the density field
        double smoothed_rho[nx][ny];

        #pragma omp parallel for collapse(2)
        for (int x = 1; x < nx - 1; x++) {
            for (int y = 1; y < ny - 1; y++) {
                smoothed_rho[x][y] = 0.25 * (lbm.rho[x][y] + lbm.rho[x+1][y] + lbm.rho[x][y+1] + lbm.rho[x+1][y+1]);
            }
        }

        // Render the smoothed density field with enhanced colour mapping
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255); // Black background
        SDL_RenderClear(renderer);

        #pragma omp parallel for collapse(2)
        for (int x = 1; x < nx - 1; x++) {
            for (int y = 1; y < ny - 1; y++) {
                // Normalise density between min_rho and max_rho
                double norm_density = (smoothed_rho[x][y] - min_rho) / (max_rho - min_rho);
                if (norm_density > 1.0) norm_density = 1.0;
                if (norm_density < 0.0) norm_density = 0.0;

                // Apply a multi-colour gradient
                int r, g, b;

                if (norm_density < 0.5) {
                    r = 0;
                    g = (int)(255 * (2 * norm_density));  // Green increases up to mid-density
                    b = (int)(255 * (1 - 2 * norm_density));  // Blue decreases to mid-density
                } else {
                    r = (int)(255 * (2 * (norm_density - 0.5)));  // Red increases after mid-density
                    g = (int)(255 * (1 - 2 * (norm_density - 0.5)));  // Green decreases after mid-density
                    b = 0;
                }

                // Rendering each square is thread-safe since each pixel is independent
                SDL_SetRenderDrawColor(renderer, r, g, b, 255);

                int cell_width = WINDOW_WIDTH / nx;
                int cell_height = WINDOW_HEIGHT / ny;

                SDL_Rect rect;
                rect.x = x * cell_width;
                rect.y = y * cell_height;
                rect.w = cell_width;
                rect.h = cell_height;

                SDL_RenderFillRect(renderer, &rect);
            }
        }

        SDL_RenderPresent(renderer);

        SDL_Delay(16);  // Delay to control the frame rate (~60 FPS)

        // Print debug information every 100 time steps
        if (t % 100 == 0) {
            printf("Time step %d: ux[0][0] = %f, uy[0][0] = %f, rho[0][0] = %f\n", t, lbm.ux[0][0], lbm.uy[0][0], lbm.rho[0][0]);
        }
    }

    // Cleanup
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    free_lbm(&lbm);
    return 0;
}
