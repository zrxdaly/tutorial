/**
 * Green2D.c - 2D Greenhouse simulation with canopy and wave-shaped roof
 * This simulation models indoor flow over vegetation canopy with a wavy roof structure
 * modified from Dai et al. (2024)
 */

#include <sys/stat.h>
#include "navier-stokes/centered.h"
#include "fractions.h"
#include "output_slices.h"
#include "view.h"

// Simulation parameters
int maxlevel, minlevel;              // Grid refinement levels
double eps;                          // Adaptation criterion
double TEND = 80.;                    // Simulation end time [s]

#include "physics.h"

// Output directory settings
char filedir[200] = "/home/dai/Documents/talks/workshop_Shanghai_oct_30/tutorial/part4_parameterization_VF/";
char W12dir[230];

// ============================================================
// GREENHOUSE ROOF PARAMETERS
// ============================================================
scalar ROOF[];                       // Roof volume fraction field (1 = inside roof, 0 = outside)
face vector fROOF[];                 // Roof face fraction field
#define NUM_WAVES 2                  // Number of waves across the domain
#define ROOF_Y_MIN 70.               // Lowest point of the roof [m]
#define ROOF_Y_MAX 100.              // Highest point of the roof [m]
#define ROOF_Y_CENTER 85.            // Center height of roof [m]
#define ROOF_AMPLITUDE 15.           // Wave amplitude [m]
#define ROOF_WAVE_NUMBER (2.*M_PI*NUM_WAVES/L0)  // Wave number [rad/m]
#define ROOF_VAL 0.5                 // Smoothing parameter for fractions

// ============================================================
// MAIN FUNCTION - Set up simulation parameters
// ============================================================
int main()
{
  // Grid configuration
  minlevel = 5;                      // Minimum refinement level
  maxlevel = 7;                      // Maximum refinement level
  N = 64;                            // Base grid size

  // Domain configuration
  L0 = 100;                          // Domain size [m]
  X0 = Y0 = 0;                       // Domain origin

  // Physics configuration
  a = av;                            // Acceleration field (buoyancy)
  mu = Km;                           // Turbulent viscosity
  Pr = unityf;                       // Prandtl number = 1
  eps = .05;                         // Adaptation criterion

  // Numerical methods
  foreach_dimension()
    u.x.refine = refine_linear;      // Momentum conserved during refinement
  p.refine = p.prolongation = refine_linear;  // Pressure interpolation
  b.gradient = minmod2;              // Flux limiter for buoyancy

  // Boundary conditions
  Boundary_C();                      // Set physics boundary conditions
  leaf_BC();                         // Set canopy boundary conditions

  run();
}

// ============================================================
// HELPER FUNCTION - Create output directory
// ============================================================
void sim_dir_create(){
    sprintf(W12dir, "%sW12/", filedir);
    if (pid() == 0)                  // Only master process creates directory
    {
        struct stat st = {0};
        if (stat(W12dir, &st) == -1)
            mkdir(W12dir, 0777);
    }
}

// ============================================================
// EVENT: INIT - Initialize flow field and variables at t=0
// ============================================================
event init(t = 0)
{
    init_physics();                  // Initialize physical fields from physics.h
    sim_dir_create();                // Create output directory

    // Add noise for initial turbulence in the canopy region
    scalar n[];
    foreach () {
        n[] = 0.;
        if (x <= 85. && x >= 15. && y < 35. && y > 65.)
            n[] = noise();
    }

    // Adapt grid based on initial conditions
    while (adapt_wavelet((scalar *){u, b, n},
           (double[]){eps, eps, 0.35 * 9.81 / 273, 0.001},
           maxlevel, minlevel).nf)
    {
        foreach () {
            b[] = STRAT(y);          // Initialize buoyancy stratification
            u.x[] = WIND(y);         // Initialize wind profile
        }
    }

    // Initialize canopy variables
    foreach() {
        TV[] = 295.15;               // Leaf temperature [K]
        e120[] = Emin;               // TKE initialization
    }
}

// ============================================================
// EVENT: ROOF - Calculate roof fraction and damp velocity
// ============================================================
event roof(i++)
{
    // Define level set for wavy roof using fractions method
    vertex scalar phi[];
    foreach_vertex() {
        // Roof height as function of x: y_roof(x) = center + amplitude * cos(k*x)
        double y_roof = ROOF_Y_CENTER + ROOF_AMPLITUDE * cos(ROOF_WAVE_NUMBER * x);

        // Level set: positive inside roof (y > y_roof), negative outside
        // Distance from point to roof surface
        phi[] = y - y_roof;
    }

    // Compute volume and face fractions from level set
    fractions(phi, ROOF, fROOF);

    // Damp velocity in roof region
    foreach() {
        if (ROOF[] > 0.5) {          // If more than 50% inside roof
            u.x[] = 0.;              // Damp horizontal velocity
            u.y[] = 0.;              // Damp vertical velocity
        }
    }

    boundary({ROOF, u});
}

// ============================================================
// EVENT: ADAPT - Adaptive mesh refinement at each time step
// ============================================================
event adapt(i++)
{
    // Add noise near canopy for better refinement
    scalar n[];
    foreach() {
        n[] = 0.;
        if (x <= 85. && x >= 15. && y < 35. && y > 65.)
            n[] = noise();
    }

    // Adapt based on velocity, buoyancy, and noise
    adapt_wavelet((scalar *){u, b, n},
                  (double[]){eps, eps, .05 * 9.81 / 273, 0.001},
                  maxlevel, minlevel);
}

// ============================================================
// EVENT: SLICE_MAKE - Output 2D slice data at end time
// ============================================================
event slice_make(t += TEND)
{
    int res = 128;
    coord slice = {1., 1., 0};
    char names[300];
    snprintf(names, 300, "%sslice_%02g", W12dir, t);
    FILE *fp = fopen(names, "w");
    output_slice(list = (scalar *){b, u.x, u.y}, fp = fp,
                 n = res, linear = true, plane = slice);
    fclose(fp);
}

// ============================================================
// EVENT: DUMP_FILE1 - Save complete simulation state for restart
// ============================================================
event dump_file1(t += TEND)
{
    char name[80];
    sprintf(name, "dump-%03d", (int)t);
    dump(file = name);
}

// ============================================================
// EVENT: OUTPUT_B - Generate visualization videos every 1 second
// ============================================================
double mov_dt = 3.;
event output_b(t += mov_dt)
{
    view(width = 1200, height = 1200, tx = -0.5, ty = -0.5);

    // Video 1: Buoyancy field with canopy and roof
    draw_vof("CL", "fCL", fc = {1,1,1});
    draw_vof("CL", "fCL", filled = 1, fc = {0,1,0});
    draw_vof("ROOF", "fROOF", lw = 2, lc = {1, 1, 1});
    draw_vof("ROOF", "fROOF", filled = 1, fc = {1, 1, 1});
    squares("b", min = 0.80, max = 0.825, map = cool_warm, linear = true);
    save("b.mp4");

    // Video 2: Horizontal velocity field with canopy and roof
    draw_vof("CL", "fCL", fc = {1,1,1});
    draw_vof("CL", "fCL", filled = 1, fc = {0,1,0});
    draw_vof("ROOF", "fROOF", lw = 2, lc = {1, 1, 1});
    draw_vof("ROOF", "fROOF", filled = 1, fc = {1, 1, 1});
    squares("u.x", min = 0, max = 0.8, map = cool_warm, linear = true);
    save("ux.mp4");

    // Video 3: Water vapor field with canopy and roof
    draw_vof("CL", "fCL", fc = {1,1,1});
    draw_vof("CL", "fCL", filled = 1, fc = {0,1,0});
    draw_vof("ROOF", "fROOF", lw = 2, lc = {1, 1, 1});
    draw_vof("ROOF", "fROOF", filled = 1, fc = {1, 1, 1});
    squares("cw", min = 0, max = 0.8, map = cool_warm, linear = true);
    save("cw.mp4");
}

// ============================================================
// EVENT: END - Simulation finished
// ============================================================
event end(t = TEND)
{
    // Simulation complete
}
