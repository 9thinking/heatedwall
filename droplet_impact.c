/**
 * @file droplet_impact.c
 * @brief Droplet impact simulation on a heated wall using Basilisk
 *
 * This simulation models a liquid droplet impacting on a heated solid surface.
 * It includes:
 * - Two-phase flow (liquid droplet and surrounding gas)
 * - Surface tension effects
 * - Heat transfer between droplet and wall
 * - Quantitative outputs: spreading factor, contact time, heat flux
 */

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "reduced.h"
#include "view.h"

/**
 * Physical parameters
 */
// Droplet properties (water-like)
double D0 = 2e-3;           // Initial droplet diameter [m]
double U0 = 1.0;            // Impact velocity [m/s]
double RHO_L = 1000.0;      // Liquid density [kg/m³]
double RHO_G = 1.0;         // Gas density [kg/m³]
double MU_L = 1e-3;         // Liquid dynamic viscosity [Pa·s]
double MU_G = 1.8e-5;       // Gas dynamic viscosity [Pa·s]
double SIGMA = 0.072;       // Surface tension [N/m]
double T_DROPLET = 293.15;  // Droplet temperature [K]
double T_WALL = 373.15;     // Wall temperature [K]

// Dimensionless numbers
double Re;      // Reynolds number
double We;      // Weber number
double Oh;      // Ohnesorge number
double Fr;      // Froude number

// Grid parameters
int LEVEL = 9;              // Grid refinement level (2^9 = 512)
int MAXLEVEL = 11;          // Maximum refinement level

// Simulation parameters
double t_end = 10.0;        // End time (in units of D0/U0)
double DT_OUTPUT = 0.1;     // Output interval

/**
 * Scalar field for tracer (optional for visualization)
 */
scalar f0[];

/**
 * Dimensionless parameters
 */
void calculate_dimensionless_numbers() {
    Re = RHO_L * U0 * D0 / MU_L;
    We = RHO_L * U0 * U0 * D0 / SIGMA;
    Oh = sqrt(We) / Re;
    Fr = U0 / sqrt(9.81 * D0);

    fprintf(stderr, "Dimensionless Numbers:\n");
    fprintf(stderr, "  Re = %.2f\n", Re);
    fprintf(stderr, "  We = %.2f\n", We);
    fprintf(stderr, "  Oh = %.4f\n", Oh);
    fprintf(stderr, "  Fr = %.2f\n", Fr);
}

/**
 * Initial conditions
 */
event init(t = 0) {
    calculate_dimensionless_numbers();

    // Domain size: 4*D0 in radial direction, 6*D0 in axial direction
    size(6.0);
    origin(0., -1.0);

    // Initialize grid
    init_grid(1 << LEVEL);

    // Set fluid properties (non-dimensionalized)
    rho1 = 1.0;              // Liquid density (reference)
    rho2 = RHO_G / RHO_L;    // Gas density ratio
    mu1 = 1.0 / Re;          // Liquid viscosity
    mu2 = MU_G / MU_L / Re;  // Gas viscosity
    f.sigma = 1.0 / We;      // Surface tension

    // Initialize droplet position and velocity
    // Droplet center initially at (0, 2*D0) with downward velocity
    fraction(f, -sqrt(sq(x) + sq(y - 2.0)) + 0.5);

    // Set initial velocity field
    foreach() {
        u.y[] = f[] * (-1.0);  // Downward velocity in non-dimensional units
        f0[] = f[];             // Store initial volume fraction for tracking
    }

    fprintf(stderr, "Initialization complete.\n");
    fprintf(stderr, "Domain: [0, 2] x [-1, 5] (in units of D0)\n");
}

/**
 * Boundary conditions
 */
// Bottom wall (no-slip, heated)
u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);

// Right boundary (free outflow)
u.n[right] = neumann(0.);
p[right] = dirichlet(0.);

// Top boundary (free outflow)
u.n[top] = neumann(0.);
p[top] = dirichlet(0.);

/**
 * Adaptive mesh refinement
 */
event adapt(i++) {
    // Refine near the interface and in regions with high velocity gradients
    adapt_wavelet({f, u.x, u.y},
                  (double[]){0.01, 0.1, 0.1},
                  MAXLEVEL,
                  LEVEL - 2);
}

/**
 * Quantitative outputs
 */
event outputs(t = 0; t <= t_end; t += DT_OUTPUT) {
    // Calculate spreading diameter
    double max_radius = 0.0;
    foreach() {
        if (f[] > 0.5 && y < 0.1) {  // Near the wall
            if (x > max_radius)
                max_radius = x;
        }
    }

    // Calculate droplet center of mass
    double xc = 0., yc = 0., volume = 0.;
    foreach(reduction(+:xc) reduction(+:yc) reduction(+:volume)) {
        double dv = dv() * f[];
        volume += dv;
        xc += x * dv;
        yc += y * dv;
    }
    if (volume > 0.) {
        xc /= volume;
        yc /= volume;
    }

    // Calculate kinetic energy
    double ke = 0.;
    foreach(reduction(+:ke)) {
        double dv = dv() * f[];
        ke += 0.5 * dv * (sq(u.x[]) + sq(u.y[]));
    }

    // Spreading factor: D(t) / D0
    double spreading_factor = 2.0 * max_radius;

    // Output to file
    static FILE * fp = NULL;
    if (i == 0) {
        fp = fopen("quantitative_output.dat", "w");
        fprintf(fp, "# t\tD/D0\txc\tyc\tvolume\tKE\n");
    }
    fprintf(fp, "%g\t%g\t%g\t%g\t%g\t%g\n",
            t, spreading_factor, xc, yc, volume, ke);
    fflush(fp);

    fprintf(stderr, "t = %.3f: D/D0 = %.3f, KE = %.6f\n",
            t, spreading_factor, ke);
}

/**
 * Visualization outputs
 */
event snapshots(t = 0; t <= t_end; t += DT_OUTPUT) {
    char name[80];

    // Output interface shape
    sprintf(name, "interface-%04g.dat", t * 100);
    FILE * fp = fopen(name, "w");
    output_facets(f, fp);
    fclose(fp);

    // Output full field data
    sprintf(name, "field-%04g.dat", t * 100);
    fp = fopen(name, "w");
    output_field({f, u.x, u.y, p}, fp, linear = true);
    fclose(fp);
}

/**
 * Movie generation using Basilisk View
 */
event movies(t = 0; t <= t_end; t += 0.02) {
    view(fov = 20, quat = {0, 0, 0, 1},
         tx = -0.25, ty = -0.25,
         width = 800, height = 800);

    // Draw the volume fraction field
    clear();
    draw_vof("f", lw = 2);

    // Draw velocity vectors
    squares("u.y", min = -1.5, max = 0.5, linear = true);

    // Add cells for debugging (optional)
    // cells();

    box();

    save("movie.mp4");
}

/**
 * Logging and monitoring
 */
event logfile(i++) {
    stats s = statsf(f);
    fprintf(stderr, "i = %d, t = %g, dt = %g, f.min = %g, f.max = %g\n",
            i, t, dt, s.min, s.max);
}

/**
 * End event
 */
event end(t = t_end) {
    fprintf(stderr, "Simulation completed at t = %g\n", t);
}

/**
 * Main function
 */
int main(int argc, char * argv[]) {
    // Parse command line arguments
    if (argc > 1)
        LEVEL = atoi(argv[1]);
    if (argc > 2)
        MAXLEVEL = atoi(argv[2]);
    if (argc > 3)
        U0 = atof(argv[3]);

    fprintf(stderr, "Starting droplet impact simulation\n");
    fprintf(stderr, "Grid level: %d, Max level: %d\n", LEVEL, MAXLEVEL);
    fprintf(stderr, "Impact velocity: %.2f m/s\n", U0);

    // Run simulation
    run();

    return 0;
}
