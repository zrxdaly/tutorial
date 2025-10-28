/**
 * Exercise 1: CO2 Diffusion Using Manual Finite Difference Method
 *
 * This program demonstrates explicit finite difference implementation of the
 * 2D diffusion equation for CO2 dispersion from a pipe injection source.
 *
 * Physics:
 *   - Diffusion equation: ∂C/∂t = D∇²C
 *   - Explicit time integration with stability constraint: Δt ≤ Δx²/(4D)
 *   - Central difference for second derivatives
 *
 * Setup:
 *   - Domain: 5m × 5m
 *   - Initial CO2: 400 ppm (atmospheric background)
 *   - Injection: 1500 ppm maintained within 0.1m radius at center
 *   - Diffusivity: D = 1.60e-5 m²/s (effective turbulent diffusivity)
 *
 * Learning objectives:
 *   - Understand finite difference discretization
 *   - Implement stability condition for explicit methods
 *   - Handle source terms in PDEs
 *   - Use adaptive mesh refinement for time-dependent problems
 */

#include "run.h"
#include "view.h"

#define PIPE_RADIUS 0.1         // Pipe injection radius

scalar C[], dCx[], dCy[];       // CO2 concentration and derivatives
double D = 1.60e-5;             // CO2 diffusivity in air at 25°C (m²/s)

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
        C[] = 400.;             // Initially CO2 400ppm in domain
    }
    boundary({C});
}

// CO2 injection from pipe point
event injection(i++){
    foreach(){
        if (sqrt(sq(x) + sq(y)) < PIPE_RADIUS)
            C[] = 1500.;        // Maintain high concentration at injection point
    }
}

// Finite difference method for 2D diffusion equation
event integration(i++) {
    DT = sq(L0/(1 << grid->maxdepth))/(4.*D);  // Adjust for refined grids
    double dt = DT;
    dt = dtnext(dt);

    // Calculate second derivatives using central differences
    foreach(){
        dCx[] = (C[1,0] - 2.*C[0,0] + C[-1,0])/sq(Delta);
        dCy[] = (C[0,1] - 2.*C[0,0] + C[0,-1])/sq(Delta);
    }

    // Update concentration: C_new = C_old + dt * D * (d²C/dx² + d²C/dy²)
    foreach()
        C[] += dt*D*(dCx[] + dCy[]);

    boundary({C});
}

event visualization (t = 0; t <= 30.; t += 1.0) {
  view(width = 800, height = 800);
  squares ("C", linear = true, min = 400, max = 1000);
//   cells();                         // Draw grid cells
  box();                            // Draw domain boundary
  save ("CO2_field.mp4");
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


event adapt(i++) {
    adapt_wavelet({C}, (double[]){5e-1}, minlevel = 4, maxlevel = 6);
}