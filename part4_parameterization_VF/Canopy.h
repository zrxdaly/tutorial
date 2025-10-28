/**
 * Canopy.h - Vegetation canopy model for greenhouse simulation
 *
 * Based on: Patton et al. (2016); Boekee et al. (2023); Dai et al. (2024)
 * This header implements:
 *   - Three cube-shaped canopy elements using fractions method
 *   - Energy balance for vegetation leaves
 *   - Convective and radiative heat exchange
 *   - Transpiration through stomatal resistance
 */

#include "fractions.h"

// ============================================================
// CANOPY GEOMETRY - Three discrete cube elements
// ============================================================
#define Zh 30.                       // Canopy top height [m]
scalar CL[];                         // Canopy volume fraction field
face vector fCL[];                   // Canopy face fraction field

#define CUBE_HEIGHT 30.              // Vertical height of each cube [m]
#define CUBE_WIDTH 5.                // Horizontal width of each cube [m]
#define CUBE_VAL 0.6                 // Smoothing parameter for fractions

// Cube center positions
#define CUBE1_X 20.                  // First cube at x=20m
#define CUBE2_X 50.                  // Second cube at x=50m
#define CUBE3_X 80.                  // Third cube at x=80m

#define PAD(s) 1.20                  // Plant area density [m^2/m^3]
#define Cd 0.20                      // Drag coefficient (Shaw & Schumann 1992)

// ============================================================
// PHYSICAL CONSTANTS
// ============================================================
#define boltz 5.67E-8                // Stefan-Boltzmann constant [W m^-2 K^-4]
#define Gconst 9.81                  // Gravitational acceleration [m s^-2]
#define T_ref 273.15                 // Reference temperature [K]
#define Kd 0.024                     // Thermal conductivity of air [W m^-1 K^-1]

// ============================================================
// RADIATION PARAMETERS
// ============================================================
#define VF_s 0.1                     // View factor of sky
#define VF_g VF_s                    // View factor of ground
#define VF_l (1. - VF_s)             // View factor of surrounding leaves
#define eps_s 0.8                    // Emissivity of sky
#define eps_g 0.98                   // Emissivity of ground
#define eps_l 0.96                   // Emissivity of leaf
#define T_s 295.15                   // Sky temperature [K]
#define T_g 295.15                   // Ground temperature [K]

// ============================================================
// LEAF PROPERTIES
// ============================================================
#define Cp_l 2.0E8                   // Leaf heat capacity [J m^-3 K^-1]
#define Cp_a 1005.                   // Air heat capacity [J kg^-1 K^-1]
#define rho_a 1.27                   // Air density at 10°C [kg m^-3]
#define dvis 1.718E-5                // Dynamic viscosity [N s m^-2]
#define vis (dvis / rho_a)           // Kinematic viscosity [m^2 s^-1]

// Leaf geometry
#define R_l 4E-1                     // Leaf radius [m]
#define L_l (2 * R_l)                // Characteristic length scale [m]
#define d_l 2.0E-4                   // Leaf thickness [m]
#define A_l (2. * M_PI * sq(R_l))    // Leaf surface area [m^2]
#define V_l (A_l / 2. * d_l)         // Leaf volume [m^3]

// ============================================================
// SCALAR FIELDS
// ============================================================
scalar TV[];                         // Leaf temperature [K]
scalar H[];                          // Convective heat flux [W m^-2]
scalar QE[];                         // Transpiration flux [kg m^-2 s^-1]
scalar Lwnet[];                      // Net longwave radiation [W m^-2]

// Transpiration parameters
double cw_sat = 1.28;                // Saturation water vapor concentration
double rs = 231.;                    // Stomatal resistance [s m^-1]

// ============================================================
// FUNCTION: leaf_BC - Set boundary conditions for canopy fields
// ============================================================
void leaf_BC()
{
    TV.refine = refine_injection;
    TV.coarsen = refine_injection;
    H.refine = refine_injection;
    H.coarsen = refine_injection;
    QE.refine = refine_injection;
    QE.coarsen = refine_injection;
    CL.refine = refine_injection;
    CL.coarsen = refine_injection;
}

// ============================================================
// EVENT: leaf_flow - Compute canopy geometry and energy balance
// ============================================================
event leaf_flow(i++)
{
    // ------------------------------------------------------------
    // STEP 1: Define canopy geometry using fractions method
    // ------------------------------------------------------------
    vertex scalar phi[];
    foreach_vertex() {
        // Define three cubes as level set functions (positive inside, negative outside)
        double cube1 = intersection(intersection(
            CUBE_HEIGHT - y,                       // Top: y < 30m
            y),                                    // Bottom: y > 0m
            CUBE_WIDTH/2. - fabs(x - CUBE1_X));   // Sides: |x-20| < 2.5m

        double cube2 = intersection(intersection(
            CUBE_HEIGHT - y,
            y),
            CUBE_WIDTH/2. - fabs(x - CUBE2_X));   // Centered at x=50m

        double cube3 = intersection(intersection(
            CUBE_HEIGHT - y,
            y),
            CUBE_WIDTH/2. - fabs(x - CUBE3_X));   // Centered at x=80m

        // Union: combine all three cubes
        phi[] = union(union(cube1, cube2), cube3);
    }

    // Compute volume and face fractions from level set
    fractions(phi, CL, fCL);
    boundary({CL});

    // ------------------------------------------------------------
    // STEP 2: Compute longwave radiation balance
    // ------------------------------------------------------------
    foreach() {
        Lwnet[] = 0.;
        if (CL[] > 0.) {
            // Incoming longwave from sky, ground, and surrounding leaves
            double Lwin = 0.5 * VF_s * eps_s * boltz * pow(T_s, 4) +
                          0.5 * VF_g * eps_g * boltz * pow(T_g, 4) +
                          1.0 * VF_l * eps_l * boltz * pow(TV[], 4);

            // Outgoing longwave from leaf surface
            double Lwout = eps_l * boltz * pow(TV[], 4);

            // Net longwave radiation
            Lwnet[] = Lwin - Lwout;
        }
    }

    // ------------------------------------------------------------
    // STEP 3: Compute convective heat exchange and transpiration
    // ------------------------------------------------------------
    foreach () {
        H[] = 0.;
        QE[] = 0.;

        if (CL[] > 0.) {
            // Air temperature from buoyancy field
            double T_a = b[] * T_ref / Gconst + T_ref;

            // Buoyancy parameter for mixed convection
            double gstar = Gconst * (TV[] - T_a) / T_a;

            // Effective velocity including buoyancy effects
            double M = sqrt(sq(u.x[]) + sq(u.y[]) + fabs(2 * L_l * gstar));

            // Reynolds and Nusselt numbers for heat transfer
            double Re = fabs(M) * L_l / vis;
            double Nu = (Re > 2E4) ? 0.032 * pow(Re, 0.8) : 0.6 * pow(Re, 0.5);

            // Aerodynamic resistance
            double rH = L_l / Nu / Kd * Cp_a * rho_a;

            // Convective heat flux [W m^-2]
            H[] = Cp_a * rho_a / rH * (TV[] - T_a);

            // Update leaf temperature (forward Euler)
            TV[] += dt * (Lwnet[] - H[]) * A_l / (Cp_l * V_l);

            // Transpiration flux [kg m^-2 s^-1]
            QE[] = (cw_sat - cw[]) / (rH + rs);
        }
    }
}