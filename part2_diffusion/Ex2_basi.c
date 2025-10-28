/**
 * Exercise 2: CO2 Diffusion Using Basilisk's Implicit Solver
 *
 * This program solves the same CO2 diffusion problem as Exercise 1, but uses
 * Basilisk's built-in implicit diffusion solver for improved stability and efficiency.
 *
 * Physics:
 *   - Diffusion equation: ∂C/∂t = D∇²C
 *   - Implicit time integration (unconditionally stable)
 *   - Built-in solver handles boundary conditions automatically
 *
 * Setup:
 *   - Same physical setup as Exercise 1
 *   - Domain: 10m × 10m
 *   - Initial CO2: 400 ppm
 *   - Injection: 1500 ppm at center (radius 0.1m)
 *   - Diffusivity: D = 0.1 m²/s
 *
 * Learning objectives:
 *   - Use Basilisk's diffusion.h solver
 *   - Understand advantages of implicit methods
 *   - Compare with explicit finite difference (Exercise 1)
 *   - Efficient handling of diffusion in complex geometries
 */

#include "run.h"
#include "diffusion.h"
#include "view.h"
#define PIPE_RADIUS 0.1  // Pipe injection radius (m)

scalar C[];  // CO2 concentration field
double D = 0.1;  // Effective turbulent diffusivity (m²/s)

int main(){
    L0 = 5.;                    // domain size: 5 m
    X0 = -L0/2.;                // Set origin at center
    Y0 = -L0/2.;
    N = 1 << 6;                 // the initial grid size is set as 2 to power of 6
    DT = sq(L0/N)/(4.*D);       // Stability condition: dt <= dx^2/(4*D)
    run();
}

event init(t = 0){
    foreach(){
        C[] = 400.;             // Background atmospheric CO2 (ppm)
    }
    boundary({C});
}

// Apply diffusion using Basilisk's implicit solver
event Diffusion(i++){
    // Adjust timestep for refined grids
    DT = sq(L0/(1 << grid->maxdepth))/(4.*D);
    double dt = DT;
    dt = dtnext(dt);
    // Solve diffusion equation: dC/dt = D * ∇²C
    const face vector kappa[] = {D, D};     // Diffusion coefficient
    diffusion(C, dt, kappa);
}

// CO2 injection source - apply after diffusion
event injection(i++){
    foreach(){
        if (sqrt(sq(x) + sq(y)) < PIPE_RADIUS){
            C[] = 1500.;        // Maintain injection concentration (ppm)
        }
    }
    boundary({C});
}

event visualization(t = 0; t <= 30.; t += 1.0){
    view(width = 800, height = 800);        // Set image size
    squares("C", linear = true, min = 400, max = 1000);  // Color map
    cells();                         // Draw grid cells
    box();
    save("CO2_field.mp4");
}

// here i add a event that can log data at the line of y = 0 using interpolate function
// do you see the output in the folder, we are going to some data processing
event printdata (t = 0; t <= 30.; t += 1.0) {
  static FILE * fp = fopen ("y0.dat","w");
  for (double x = -L0/2; x<L0/2; x+=L0/N){
    fprintf (fp, "%g %g %g\n", t, x, interpolate(C, x, 0));}
  fprintf (fp, "\n");
  fflush (fp);
}

// Adaptive mesh refinement
event adapt(i++){
    adapt_wavelet({C}, (double[]){5e-1}, minlevel = 4, maxlevel = 6);
}