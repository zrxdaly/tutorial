/**
 * Bonus Challenge 2: Circular Leaf Shape
 *
 * This program compares CO2 flux between circular and elliptical leaves.
 *
 * Challenge:
 *   - Replace ellipse with circle of same surface area
 *   - Compare flux between circle and ellipse
 *   - Understand shape effects on mass transfer
 *
 * Key Questions:
 *   - Is the ellipse more efficient than a circle?
 *   - Does shape matter more with or without wind?
 *   - What is the optimal leaf shape for maximum CO2 uptake?
 *
 * Note:
 *   - Ellipse area = π * r1 * r2 = π * 5 * 1 ≈ 15.7
 *   - Circle area = π * r² → r = √5 ≈ 2.24 (for same area)
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
double wind_in = 1.;        // With wind to see shape effect
double s_in = 40.0;
double s_ls = 20.0;

// ============================================================================
// BOUNDARY CONDITIONS
// ============================================================================
u.t[embed] = dirichlet (0.);
u.n[embed] = dirichlet (0.);
s[embed]   = dirichlet (s_ls);

// ============================================================================
// LEAF GEOMETRY - CIRCULAR LEAF
// ============================================================================
// Circle with same area as ellipse (5 × 1)
double r_circle = 2.236;    // √5 ≈ 2.236 for equivalent area

#define CIRCLE (sq(x) + sq(y) - sq(r_circle))

// For comparison, original ellipse would be:
// double r1 = 5, r2 = 1;
// #define ELLIPSE (sq(x/r1) + sq(y/r2) - 1.)

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

// ============================================================================
// INITIALIZATION - Use CIRCLE instead of ELLIPSE
// ============================================================================
event init (t = 0) {
  // Refine mesh around circular leaf
  refine (CIRCLE <  2.5 && level  <  maxlevel - 1);
  refine (CIRCLE > -0.5 && CIRCLE <  0.5 && level  <  maxlevel);

  vertex scalar phi[];
  foreach_vertex()
    phi[] = CIRCLE;    // Use circular geometry
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
