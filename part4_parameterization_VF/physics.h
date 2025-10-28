/**
 * physics.h - Physical parameterizations for indoor climate simulation
 *
 * This header implements:
 *   - Physical constants and initial conditions
 *   - Boundary conditions for velocity, buoyancy, and tracers
 *   - Buoyancy-driven acceleration (gravity forcing)
 *   - Canopy drag effects on momentum
 *   - Scalar diffusion with source terms from vegetation
 */

// ============================================================
// PHYSICAL CONSTANTS
// ============================================================
#define CP 1005.                     // Specific heat of air at constant pressure [J kg^-1 K^-1]
#define gCONST 9.81                  // Gravitational acceleration [m s^-2]
#define TREF 273.15                  // Reference temperature [K]
#define karman 0.4                   // von Kármán constant (dimensionless)

// Surface roughness parameters
#define roughY0u 0.1                 // Roughness length for momentum [m]
#define roughY0h 0.1                 // Roughness length for heat [m]

// ============================================================
// INITIAL CONDITIONS
// ============================================================
#define U0 0.2                       // Reference wind speed [m s^-1]
#define WIND(s) (0.5)                // Initial wind profile [m s^-1]
#define BSURF (gCONST / TREF * 22.5) // Surface buoyancy [m s^-2]
#define STRAT(s) (gCONST / TREF * 22.5) // Buoyancy stratification [m s^-2]

// ============================================================
// SCALAR FIELDS AND TRACERS
// ============================================================
scalar b[];                          // Buoyancy field [m s^-2]
scalar U[];                          // Velocity magnitude [m s^-1]
scalar e120[];                       // Turbulent kinetic energy (TKE) [m^2 s^-2]
scalar cw[];                         // Water vapor concentration [kg kg^-1]
scalar *tracers = {b, e120, cw};     // List of advected tracers

// ============================================================
// CANOPY AND TURBULENCE MODELS
// ============================================================
#define CANOPY 1                     // Enable canopy model (1=on, 0=off)
#if CANOPY
#include "Canopy.h"
#endif

#include "SGS_TKE.h"                 // Sub-grid scale turbulence model

// ============================================================
// GLOBAL VECTORS
// ============================================================
face vector av[];                    // Acceleration vector (buoyancy + drag) [m s^-2] 

// ============================================================
// FUNCTION: Boundary_C - Set boundary conditions for all fields
// ============================================================
void Boundary_C(){
    b.nodump = false;                // Allow buoyancy to be saved in dumps

    // ------------------------------------------------------------
    // Velocity boundary conditions
    // ------------------------------------------------------------
    // Bottom wall: no-slip condition
    u.n[bottom] = dirichlet(0.);     // Normal velocity = 0
    u.t[bottom] = dirichlet(0.);     // Tangential velocity = 0

    // Top wall: no-slip condition
    u.n[top] = dirichlet(0.);
    u.t[top] = dirichlet(0.);

    // Left/right: periodic (flow wraps around)
    periodic(left);

    // ------------------------------------------------------------
    // Buoyancy boundary conditions
    // ------------------------------------------------------------
    b[bottom] = dirichlet(STRAT(y)); // Fixed buoyancy at bottom
    b[top] = dirichlet(STRAT(y));    // Fixed buoyancy at top

    // ------------------------------------------------------------
    // TKE boundary conditions
    // ------------------------------------------------------------
    Evis[bottom] = dirichlet(0.);    // Zero TKE at walls (resolved boundary layer)
    Evis[top] = dirichlet(0.);
}

// ============================================================
// FUNCTION: init_physics - Initialize flow field variables
// ============================================================
void init_physics(){
    foreach () {
        cw[] = 0.6;                  // Initial water vapor concentration [kg kg^-1]
        b[] = STRAT(y);              // Initial buoyancy stratification [m s^-2]
        u.x[] = WIND(y);             // Initial horizontal velocity [m s^-1]
        u.y[] = 0.;                  // Initial vertical velocity [m s^-1]
    }
}

// ============================================================
// EVENT: acceleration - Compute acceleration from buoyancy and canopy drag
// ============================================================
event acceleration(i++)
{
    // ------------------------------------------------------------
    // STEP 1: Buoyancy forcing (gravity)
    // ------------------------------------------------------------
    foreach_face(y) {
        // Vertical acceleration from buoyancy (averaged to face)
        av.y[] = (b[] + b[0,-1])/2.;
    }

    // ------------------------------------------------------------
    // STEP 2: Canopy drag (momentum sink)
    // ------------------------------------------------------------
    #if CANOPY
    // Calculate velocity magnitude at cell centers
    foreach () {
        U[] = sqrt(sq(u.x[]) + sq(u.y[]));
    }

    // Apply drag force to horizontal faces
    // Drag = -Cd * PAD * |U| * u * CL
    foreach_face() {
        // Interpolate velocity magnitude and velocity to face
        double U_face = (U[] + U[-1])/2.;
        double u_face = (u.x[] + u.x[-1])/2.;
        double CL_face = (CL[] + CL[-1])/2.;

        // Add drag acceleration (negative = opposing flow)
        av.x[] = av.x[] - Cd * PAD(y) * U_face * u_face * CL_face;
    }
    #endif
}

// ============================================================
// EVENT: tracer_diffusion - Solve diffusion with source terms from canopy
// ============================================================
mgstats mgb;                         // Multigrid solver statistics

event tracer_diffusion(i++)
{
    scalar r[];                      // Source term for buoyancy [m s^-3]
    scalar r_cw[];                   // Source term for water vapor [kg kg^-1 s^-1]

    foreach() {
        r[] = 0;
        r_cw[] = 0;

        #if CANOPY
        // Source from sensible heat flux (converts to buoyancy)
        // H [W m^-2] -> [K m s^-1] -> [m s^-3] (buoyancy units)
        r[] = r[] + H[] / (Cp_a * rho_a) * (gCONST / TREF) * PAD(y) * CL[];

        // Source from transpiration (converts to concentration change)
        // QE [kg m^-2 s^-1] -> [kg kg^-1 s^-1] (concentration units)
        r_cw[] = r_cw[] + QE[] * PAD(y) * CL[];
        #endif
    }

    // Alternative: surface flux calculations (commented out)
    // double flx = 0, bt = 0;
    // double fctr = rho_a * CP * TREF / gCONST;
    // foreach_boundary(bottom, reduction(+:flx) reduction(+:bt)) {
    //     flx = flx + (QFLX + GFLX) * sq(Delta);
    //     bt = bt + BSURF * sq(Delta);
    // }
    // bt = bt/sq(L0);
    // flx = flx/sq(L0);
    // fprintf(stderr, "soil=%g %g %g %d\n", t, fctr*flx, bt * TREF / gCONST, i);

    // Solve diffusion equation with source terms
    mgb = diffusion(b, dt, Kh, r = r);      // Buoyancy diffusion
    diffusion(cw, dt, Kh, r = r_cw);        // Water vapor diffusion
}
