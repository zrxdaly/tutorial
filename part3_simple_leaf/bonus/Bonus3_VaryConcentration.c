/**
 * Bonus Challenge 3: Vary Concentration Gradient
 *
 * This program investigates whether CO2 flux is linear with concentration difference.
 *
 * Challenge:
 *   - Test different ambient (s_in) and leaf surface (s_ls) concentrations
 *   - Measure flux vs. concentration gradient
 *   - Verify Fick's law: Flux ∝ ΔC
 *
 * Key Questions:
 *   - Is flux linear with concentration difference?
 *   - What happens with very large or very small gradients?
 *   - Does linearity depend on Re or leaf orientation?
 *
 * Test Cases:
 *   Case 1 (default): ΔC = 40 - 20 = 20 mmol/m³
 *   Case 2 (high):    ΔC = 60 - 10 = 50 mmol/m³
 *   Case 3 (low):     ΔC = 30 - 25 = 5 mmol/m³
 */

#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"

scalar s[], * tracers = {s};

// ============================================================================
// SIMULATION PARAMETERS
// ============================================================================
int maxlevel = 9;
double ue = 0.05;
double be = 0.05;
double wind_in = 1.;

// ============================================================================
// BONUS CHALLENGE: Try different concentration gradients
// ============================================================================
// Uncomment one case to test:

// Case 1: Default gradient (ΔC = 20)
double s_in = 40.0;         // Ambient CO2
double s_ls = 20.0;         // Leaf surface CO2

// Case 2: Large gradient (ΔC = 50) - Uncomment to test
// double s_in = 60.0;         // Higher ambient CO2
// double s_ls = 10.0;         // Lower leaf surface CO2

// Case 3: Small gradient (ΔC = 5) - Uncomment to test
// double s_in = 30.0;         // Lower ambient CO2
// double s_ls = 25.0;         // Higher leaf surface CO2

// Case 4: Very large gradient (ΔC = 80) - Uncomment to test
// double s_in = 100.0;        // Very high ambient (CO2 enrichment)
// double s_ls = 20.0;         // Normal leaf surface

// ============================================================================
// BOUNDARY CONDITIONS
// ============================================================================
u.t[embed] = dirichlet (0.);
u.n[embed] = dirichlet (0.);
s[embed]   = dirichlet (s_ls);

// ============================================================================
// LEAF GEOMETRY
// ============================================================================
double r1 = 5, r2 = 1;      // Horizontal ellipse
#define ELLIPSE (sq(x/r1) + sq(y/r2) - 1.)

// ============================================================================
face vector muc[];
double Re = 50;

int main() {
  periodic (left);
  L0 = 120;
  X0 = Y0 = -L0/2;
  mu = muc;
  run();
}

event properties (i++) {
  foreach_face()
    muc.x[] = fs.x[]/Re;
  boundary ((scalar*){muc});
}

event init (t = 0) {
  refine (ELLIPSE <  2.5 && level  <  maxlevel - 1);
  refine (ELLIPSE > -0.5 && ELLIPSE <  0.5 && level  <  maxlevel);
  vertex scalar phi[];
  foreach_vertex()
    phi[] = ELLIPSE;
  boundary ({phi});
  fractions (phi, cs, fs);
  foreach(){
    s[] = cs[] > 0? s_in: 0.;
    u.x[] = cs[] > 0? wind_in: 0.;
  }
  boundary ({s, u.x});
}

event tracer_diffusion (i++) {
  diffusion (s, dt, muc);
}

event force (i++) {
  double FB = L0/5., tau = 1;
  foreach() {
    if (x < X0 + FB) {
      s[] -= (s[] - s_in)*dt/tau;
      u.y[] -= u.y[]*dt/tau;
      u.x[] -= (u.x[] - wind_in)*dt/tau;
    }
  }
  boundary ({s, u});
}

event adapt (i++) {
  adapt_wavelet ({cs, s, u}, (double[]){5e-1, be, ue,ue}, maxlevel, 5);
}

event mov (t += 1) {
  scalar m[];
  foreach()
    m[] = cs[] - 0.5;
  boundary ({m});
  // Adjust color scale based on concentration range
  output_ppm (s, file = "s.mp4", n = 512, mask = m,
	      linear = true, max = s_in, min = s_ls, map = cool_warm,
	      box = {{X0 + 15, -15},{X0 + L0, 15}});
}

event diag_flux (t = 300; t += 1) {
  double flx = 0;
  foreach(reduction(+:flx)) {
    double val = 0, e = embed_flux (point, s, mu, &val);
    if (val)
      flx += (val - e*s[])*sq(Delta);
  }
  static FILE * fp = fopen ("diag1", "w");
  // Output: time, flux, concentration_difference
  fprintf (fp, "%g %g %g\n", t, flx, s_in - s_ls);
}

event stop (t = 400);
