/**
 * Bonus Challenge 4: Inclined Leaf
 *
 * This program simulates CO2 flux to a leaf inclined at various angles.
 *
 * Challenge:
 *   - Rotate the ellipse by different angles (0°, 45°, 90°)
 *   - Find the optimal inclination angle for maximum flux
 *   - Understand angle-wind interaction
 *
 * Key Questions:
 *   - What angle gives maximum CO2 flux?
 *   - Is 45° always optimal?
 *   - How does optimal angle depend on Re and aspect ratio?
 *
 * Implementation:
 *   - Use coordinate rotation: (x', y') = (x*cos(θ) + y*sin(θ), -x*sin(θ) + y*cos(θ))
 *   - Apply rotation to ellipse equation
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
double wind_in = 1.;        // Horizontal wind
double s_in = 40.0;
double s_ls = 20.0;

// ============================================================================
// BOUNDARY CONDITIONS
// ============================================================================
u.t[embed] = dirichlet (0.);
u.n[embed] = dirichlet (0.);
s[embed]   = dirichlet (s_ls);

// ============================================================================
// LEAF GEOMETRY - INCLINED ELLIPSE
// ============================================================================
double r1 = 5, r2 = 1;      // Ellipse semi-axes

// ============================================================================
// BONUS CHALLENGE: Try different inclination angles
// ============================================================================
// Uncomment one angle to test:

double theta = 0.0;         // 0° - Horizontal (same as Ex2)
// double theta = M_PI/6;      // 30° inclination
// double theta = M_PI/4;      // 45° inclination
// double theta = M_PI/3;      // 60° inclination
// double theta = M_PI/2;      // 90° - Vertical (same as Ex4)

// Rotated ellipse equation:
// Original: (x/r1)² + (y/r2)² = 1
// Rotated:  ((x*cos + y*sin)/r1)² + ((-x*sin + y*cos)/r2)² = 1
#define ELLIPSE_ROTATED (sq((x*cos(theta) + y*sin(theta))/r1) + \
                         sq((-x*sin(theta) + y*cos(theta))/r2) - 1.)

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
  // Refine around rotated ellipse
  refine (ELLIPSE_ROTATED <  2.5 && level  <  maxlevel - 1);
  refine (ELLIPSE_ROTATED > -0.5 && ELLIPSE_ROTATED <  0.5 && level  <  maxlevel);

  vertex scalar phi[];
  foreach_vertex()
    phi[] = ELLIPSE_ROTATED;    // Use rotated ellipse
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
  // Output: time, flux, angle_in_degrees
  fprintf (fp, "%g %g %g\n", t, flx, theta*180./M_PI);
}

event stop (t = 400);
