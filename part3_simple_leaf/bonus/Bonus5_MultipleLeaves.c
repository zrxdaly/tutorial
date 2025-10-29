/**
 * Bonus Challenge 5: Multiple Leaves
 *
 * This program simulates CO2 flux to multiple leaves to study interaction effects.
 *
 * Challenge:
 *   - Add a second leaf downstream
 *   - Study shading and wake effects
 *   - Compare total flux to sum of individual leaves
 *
 * Key Questions:
 *   - Does the downstream leaf get less CO2 due to upstream depletion?
 *   - What is the optimal spacing between leaves?
 *   - Is total flux sub-additive due to interaction?
 *
 * Configuration:
 *   - Leaf 1 (upstream): centered at x = -10
 *   - Leaf 2 (downstream): centered at x = +10
 *   - Both horizontal ellipses with r1=5, r2=1
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
double wind_in = 1.;        // Horizontal wind to see interaction
double s_in = 40.0;
double s_ls = 20.0;

// ============================================================================
// BOUNDARY CONDITIONS
// ============================================================================
u.t[embed] = dirichlet (0.);
u.n[embed] = dirichlet (0.);
s[embed]   = dirichlet (s_ls);

// ============================================================================
// LEAF GEOMETRY - TWO LEAVES
// ============================================================================
double r1 = 5, r2 = 1;      // Ellipse semi-axes

// ============================================================================
// BONUS CHALLENGE: Adjust leaf positions and spacing
// ============================================================================
// Leaf 1 (upstream) position
double X1 = -10.0;          // x-coordinate of first leaf center
double Y1 = 0.0;            // y-coordinate of first leaf center

// Leaf 2 (downstream) position
double X2 = 10.0;           // x-coordinate of second leaf center (spacing = 20)
double Y2 = 0.0;            // y-coordinate of second leaf center

// Try different spacings:
// double X2 = 8.0;   // Close spacing (spacing = 18)
// double X2 = 15.0;  // Wide spacing (spacing = 25)
// double Y2 = 5.0;   // Offset vertically

// Define two ellipses
#define LEAF1 (sq((x-X1)/r1) + sq((y-Y1)/r2) - 1.)
#define LEAF2 (sq((x-X2)/r1) + sq((y-Y2)/r2) - 1.)

// Combined geometry: solid where either leaf is present
#define BOTH_LEAVES (min(LEAF1, LEAF2))

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
  // Refine around both leaves
  refine (BOTH_LEAVES <  2.5 && level  <  maxlevel - 1);
  refine (BOTH_LEAVES > -0.5 && BOTH_LEAVES <  0.5 && level  <  maxlevel);

  vertex scalar phi[];
  foreach_vertex()
    phi[] = BOTH_LEAVES;    // Use combined geometry
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

// Calculate flux for each leaf separately
event diag_flux (t = 300; t += 1) {
  double flx_total = 0;     // Total flux (both leaves)
  double flx_leaf1 = 0;     // Flux to leaf 1 only
  double flx_leaf2 = 0;     // Flux to leaf 2 only

  foreach(reduction(+:flx_total) reduction(+:flx_leaf1) reduction(+:flx_leaf2)) {
    double val = 0, e = embed_flux (point, s, mu, &val);
    if (val) {
      double flux_contribution = (val - e*s[])*sq(Delta);
      flx_total += flux_contribution;

      // Determine which leaf this flux belongs to
      // Check if point is closer to leaf1 or leaf2
      double dist1 = sqrt(sq(x - X1) + sq(y - Y1));
      double dist2 = sqrt(sq(x - X2) + sq(y - Y2));

      if (dist1 < dist2) {
        flx_leaf1 += flux_contribution;
      } else {
        flx_leaf2 += flux_contribution;
      }
    }
  }

  static FILE * fp = fopen ("diag1", "w");
  // Output: time, total_flux, leaf1_flux, leaf2_flux
  fprintf (fp, "%g %g %g %g\n", t, flx_total, flx_leaf1, flx_leaf2);
}

event stop (t = 400);
