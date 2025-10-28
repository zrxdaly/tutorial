/**
 * Bonus Challenge 1: Vary Reynolds Number
 *
 * This program investigates how Reynolds number affects CO2 flux and the
 * interaction between flow regime and leaf orientation.
 *
 * Challenge:
 *   - Modify Re to 10, 50, 100, 500
 *   - Compare flux-orientation relationship at different Re
 *   - Understand transition from viscous to inertial flow regimes
 *
 * Key Questions:
 *   - Does higher Re always mean higher flux?
 *   - How does Re affect the boundary layer thickness?
 *   - At what Re does orientation matter most?
 */

#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"

scalar s[], * tracers = {s};

// ============================================================================
// SIMULATION PARAMETERS - MODIFY Re HERE
// ============================================================================
int maxlevel = 9;
double ue = 0.05;
double be = 0.05;
double wind_in = 1.;        // With wind to see Re effect
double s_in = 40.0;
double s_ls = 20.0;

// ============================================================================
// BONUS CHALLENGE: Try different Reynolds numbers
// ============================================================================
// Uncomment one of the following to test different flow regimes:

double Re = 10;      // Very viscous, thick boundary layer
// double Re = 50;      // Moderate (default case)
// double Re = 100;     // Transitional
// double Re = 500;     // More inertial, thin boundary layer

// ============================================================================
// BOUNDARY CONDITIONS
// ============================================================================
u.t[embed] = dirichlet (0.);
u.n[embed] = dirichlet (0.);
s[embed]   = dirichlet (s_ls);

// ============================================================================
// LEAF GEOMETRY - Try both horizontal and vertical
// ============================================================================
// Horizontal leaf
double r1 = 5, r2 = 1;

// Vertical leaf (uncomment to test)
// double r1 = 1, r2 = 5;

#define ELLIPSE (sq(x/r1) + sq(y/r2) - 1.)

// ============================================================================
face vector muc[];

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
  output_ppm (s, file = "s.mp4", n = 512, mask = m,
	      linear = true, max = 40., min = 20., map = cool_warm,
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
  fprintf (fp, "%g %g\n", t, flx);
}

event stop (t = 400);
