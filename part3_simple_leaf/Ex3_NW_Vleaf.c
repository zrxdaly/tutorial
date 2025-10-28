/**
 * Exercise 1: CO2 Scalar Transfer with an Elliptical Leaf
 *
 * This program simulates CO2 diffusive flux from the outer environment to a leaf
 * using Navier-Stokes flow solver with embedded boundary method and scalar transport.
 *
 * Physics:
 *   - Navier-Stokes equations for incompressible flow
 *   - Scalar transport equation: ∂s/∂t + u·∇s = ∇·(D∇s)
 *   - Embedded boundary method for the leaf geometry
 *   - No-slip boundary condition at leaf surface
 *
 * Setup:
 *   - Domain: 120 × 120 units, centered at origin
 *   - Leaf geometry: Ellipse with r_major = 5, r_minor = 1
 *   - Flow: Re = 50, horizontal wind = 1 m/s
 *   - CO2 concentration: ambient = 40 mmol/m³, leaf surface = 20 mmol/m³
 *   - Concentration gradient drives diffusive flux into leaf
 *
 * Learning objectives:
 *   - Effect of wind & leaf orientation on gas flux
 *   - Understand coupled flow and scalar transport
 *   - Learn embedded boundary method for complex geometries
 */
// ============================================================================
// INCLUDE LIBRARIES
// ============================================================================
// These libraries provide the necessary functions for our simulation
#include "embed.h"                      // For embedded boundary (the leaf)
#include "navier-stokes/centered.h"     // For fluid flow equations
#include "tracer.h"                     // For tracking scalar transport (CO2)
#include "diffusion.h"                  // For diffusion of CO2

// ============================================================================
// SCALAR FIELD DEFINITION
// ============================================================================
// Define the CO2 concentration field 's' and register it as a tracer
scalar s[], * tracers = {s};

// ============================================================================
// SIMULATION PARAMETERS
// ============================================================================
int maxlevel = 9;           // Maximum grid refinement level (higher = finer mesh)
double ue = 0.05;           // Error tolerance for velocity adaptation
double be = 0.05;           // Error tolerance for scalar adaptation
double wind_in = 0.;        // Incoming wind velocity (horizontal) m s-1
double s_in = 40.0;         // average CO2 concentration in the incoming air (high) mmol m-3
double s_ls = 20.0;         // average CO2 concentration at the leaf surface (low)  mmol m-3

// ============================================================================
// BOUNDARY CONDITIONS ON THE LEAF SURFACE
// ============================================================================
// The leaf is a no-slip boundary (velocity = 0) with fixed CO2 concentration
u.t[embed] = dirichlet (0.);      // Tangential velocity = 0 (no-slip)
u.n[embed] = dirichlet (0.);      // Normal velocity = 0 (no penetration)
s[embed]   = dirichlet (s_ls);    // CO2 concentration at leaf surface = s_ls

// ============================================================================
// LEAF GEOMETRY
// ============================================================================
// reverse the major and minor axis 
double r1 = 1, r2 = 5;
#define ELLIPSE (sq(x/r1) + sq(y/r2) - 1.)  // function describes the leaf ellipse (0 at surface)

// ============================================================================
// FLUID PROPERTIES
// ============================================================================
face vector muc[];          // Kinematic viscosity field
double Re = 50;             // Reynolds number (Re = U*L/nu)

// ============================================================================
// MAIN FUNCTION - SETUP THE SIMULATION DOMAIN
// ============================================================================
int main() {
  periodic (left);          // Periodic boundary on left side (flow wraps around)
  L0 = 120;                 // Domain size: 120 x 120 units
  X0 = Y0 = -L0/2;          // Center the domain at origin (leaf at center)
  mu = muc;                 // Assign viscosity field
  run();                    // Start the simulation
}

// ============================================================================
// EVENT: UPDATE FLUID PROPERTIES
// ============================================================================
// This event runs every time step to update the kinematic viscosity
// The viscosity is scaled by the volume fraction to account for the embedded leaf
// no worries about this part
event properties (i++) {
  foreach_face()
    muc.x[] = fs.x[]/Re;      // Viscosity = volume_fraction / Reynolds_number
  boundary ((scalar*){muc});   // Update boundary conditions for viscosity
}

// this event is the initial setup and condition for the simulation
event init (t = 0) {
  // refine the local grid
  refine (ELLIPSE <  2.5 && level  <  maxlevel - 1);
  refine (ELLIPSE > -0.5 && ELLIPSE <  0.5 && level  <  maxlevel);
  // initially calculate the elliptical leaf shape
  vertex scalar phi[];
  foreach_vertex()
    phi[] = ELLIPSE;
  boundary ({phi});
  fractions (phi, cs, fs);
  // the initial condition for CO2 concentration and velocity
  foreach(){
    s[] = cs[] > 0? s_in: 0.;
    u.x[] = cs[] > 0? wind_in: 0.;
  }
  boundary ({s, u.x});
}

// the event for CO2 diffusion solver
event tracer_diffusion (i++) {
  diffusion (s, dt, muc);
}

// don't worry about this part, it is just a trick to fix the boundary condition
// for the CO2 concentration and velocity
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

// this is the event for adaptive mesh refinement
// it will refine the grid based on the error tolerance
event adapt (i++) {
  adapt_wavelet ({cs, s, u}, (double[]){5e-1, be, ue,ue}, maxlevel, 5);
}

// Generate video frames of CO2 concentration field every time unit
event mov (t += 1) {
  scalar m[];                         // Mask to hide leaf interior in visualization
  foreach()
    m[] = cs[] - 0.5;                 // Positive in fluid, negative in leaf
  boundary ({m});
  output_ppm (s, file = "s.mp4", n = 512, mask = m,
	      linear = true, max = 40., min = 20., map = cool_warm,
	      box = {{X0 + 15, -15},{X0 + L0, 15}});
  // Creates video showing CO2 from blue (20) to red (40) mmol/m3
}

// This event calculates the theoretical CO2 flux through the leaf surface
// This is the KEY OUTPUT: how much CO2 the leaf maximally can absorb from the air
// Don't worry about the details of how it is calculated
event diag_flux (t = 300; t += 1) {
  double flx = 0;                       // Total CO2 flux accumulator
  foreach(reduction(+:flx)) {           // Loop and sum over all grid cells
    double val = 0, e = embed_flux (point, s, mu, &val);
    if (val)                            // If cell contains part of leaf boundary
      flx += (val - e*s[])*sq(Delta);   // Add CO2 flux at this location
  }
  static FILE * fp = fopen ("diag1", "w");
  fprintf (fp, "%g %g\n", t, flx);    // Write: time, total_CO2_flux
}

event stop (t = 400);                   // End simulation at t = 400