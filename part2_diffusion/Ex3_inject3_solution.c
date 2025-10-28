/**
 * Exercise 3: Multiple CO2 Injection Sources
 *
 * This program extends Exercise 2 by introducing three CO2 injection nozzles
 * to study the interaction and superposition of multiple diffusion sources.
 *
 * Physics:
 *   - Same diffusion equation: ∂C/∂t = D∇²C
 *   - Multiple source terms (three injection points)
 *   - Plume interaction and merging dynamics
 *
 * Setup:
 *   - Three nozzles at: (0, 0), (-1.25, 0), (1.25, 0)
 *   - Each maintains 1500 ppm within 0.1 m radius
 *   - Background: 400 ppm
 *   - Observe how plumes spread and merge
 *
 * New features:
 *   - Data export to file for quantitative analysis
 *   - Centerline concentration profiles over time
 *   - Can be used to validate numerical solutions
 *
 * Learning objectives:
 *   - Handle multiple source terms
 *   - Export simulation data for post-processing
 *   - Analyze spatial and temporal evolution
 *   - Prepare for Exercise 4 (optimization task)
 */

#include "run.h"
#include "diffusion.h"
#include "view.h"

scalar C[];  // CO2 concentration field
double D = 0.1;  // Effective turbulent diffusivity (m²/s)
#define PIPE_RADIUS 0.1  // Pipe injection radius (m)

int main(){
    L0 = 5.;
    X0 = -L0/2.;
    Y0 = -L0/2.;
    N = 1 << 6;
    DT = sq(L0/N)/(4.*D);  // Initial timestep based on stability
    run();
}

event init(t = 0){
    foreach(){
        C[] = 400.;  // Background atmospheric CO2 (ppm)
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
    const face vector kappa[] = {D, D};  // Diffusion coefficient
    diffusion(C, dt, kappa);
}

// CO2 injection source - apply after diffusion
event injection(i++){
    double NOZZLE1_X = 0.0;
    double NOZZLE1_Y = 0.0;
    double NOZZLE2_X = -1.25;
    double NOZZLE2_Y = 0;
    double NOZZLE3_X = 1.25;
    double NOZZLE3_Y = 0;
    foreach(){
        // Nozzle 1 (center)
        double dist1 = sqrt(sq(x - NOZZLE1_X) + sq(y - NOZZLE1_Y));
        // Nozzle 2 (left side)
        double dist2 = sqrt(sq(x - NOZZLE2_X) + sq(y - NOZZLE2_Y));
        // Nozzle 3 (right side)
        double dist3 = sqrt(sq(x - NOZZLE3_X) + sq(y - NOZZLE3_Y));
        if (dist1 < PIPE_RADIUS || dist2 < PIPE_RADIUS || dist3 < PIPE_RADIUS){
            C[] = 1500.;  // Maintain injection concentration (ppm)
        }
    }
    boundary({C});
}

event visualization(t = 0; t <= 10.; t += 0.3){
    view(width = 800, height = 800);
    squares("C", linear = true, min = 400, max = 1000);
    box();
    save("CO2_field.mp4");
}

event printdata (t = 0; t <= 10.; t += 0.3) {
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