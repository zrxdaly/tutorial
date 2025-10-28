/**
 * SGS_TKE.h - Sub-grid scale turbulence model based on TKE equation
 *
 * This header implements a 1.5-order turbulence closure:
 *   - Prognostic equation for turbulent kinetic energy (TKE)
 *   - Diagnostic mixing length based on stability
 *   - Eddy viscosity and diffusivity from TKE and mixing length
 *   - Shear production, buoyancy, dissipation, and canopy drag effects
 *
 * Based on Deardorff (1980) and modifications for canopy flows (Dai et al., 2012) 
 */

#include "tracer.h"
#include "diffusion.h"

// ============================================================
// SGS TURBULENCE PARAMETERS
// ============================================================
#define Emin 5.0E-5                  // Minimum TKE (numerical stability) [m^2 s^-2]
#define kappa 0.4                    // von Kármán constant (dimensionless)
#define ce1 0.19                     // Dissipation coefficient 1 (dimensionless)
#define ce2 0.51                     // Dissipation coefficient 2 (dimensionless)
#define Cm 0.12                      // Momentum diffusion coefficient (dimensionless)
#define Cn 0.76                      // Stability coefficient (dimensionless)
#define eps1 1.0E-10                 // Minimum buoyancy gradient (numerical stability) [s^-2]

// ============================================================
// TURBULENT DIFFUSIVITIES AND FIELDS
// ============================================================
face vector Km[];                    // Eddy viscosity (momentum) [m^2 s^-1]
face vector Kh[];                    // Eddy diffusivity (heat/scalars) [m^2 s^-1]
face vector Ke[];                    // Eddy diffusivity (TKE) [m^2 s^-1]
(const) face vector Pr;              // Turbulent Prandtl number (dimensionless)
scalar Evis[];                       // Cell-centered eddy viscosity [m^2 s^-1]
scalar e12p[];                       // TKE source/sink term [m^2 s^-3]
// ============================================================
// EVENT: Km_tke - Compute eddy viscosity and TKE source terms
// ============================================================
event Km_tke(i++)
{
  foreach ()
  {
    // ------------------------------------------------------------
    // STEP 1: Calculate mixing length based on stability
    // ------------------------------------------------------------
    double lambda;

    // Vertical buoyancy gradient (central difference)
    double dbdz = (b[0, 1] - b[0, -1]) / (Delta * 2);

    // Limit gradient for numerical stability
    if (dbdz < 0 && dbdz >= -eps1)
      dbdz = -eps1;
    if (dbdz >= 0 && dbdz <= eps1)
      dbdz = eps1;

    // Mixing length depends on stability
    if (dbdz <= 0) {
      // Unstable: mixing length = grid size (neutral case)
      lambda = Delta;
    }
    else {
      // Stable: mixing length limited by buoyancy
      // λ = (κy^-1 + (Cn*e/√N²)^-1)^-1
      e120[] = (e120[] <= Emin) ? Emin : e120[];
      lambda = pow(pow(kappa * y, -1) + pow(Cn * e120[] / sqrt(fabs(dbdz)), -1), -1);
    }

    // Eddy viscosity: Km = Cm * λ * √e
    Evis[] = Cm * lambda * e120[];

    // ------------------------------------------------------------
    // STEP 2: Calculate strain rate tensor (S_ij * S_ij)
    // ------------------------------------------------------------
    // Diagonal components (normal strain rates)
    double tdef2 = 2 / sq(Delta) * (sq(u.x[1, 0, 0] - u.x[]) +
                                     sq(u.y[0, 1, 0] - u.y[]) +
                                     sq(u.z[0, 0, 1] - u.z[]));

    // Off-diagonal components (shear strain rates)
    // xy components
    tdef2 += 0.25 * (sq((u.y[0, 0, 1] - u.y[-1, 0, 1]) / Delta + (u.x[0, 0, 1] - u.x[0, 0, 0]) / Delta) +
                    sq((u.y[0, 0, 0] - u.y[-1, 0, 0]) / Delta + (u.x[0, 0, 0] - u.x[0, 0, -1]) / Delta) +
                    sq((u.y[1, 0, 0] - u.y[0, 0, 0]) / Delta + (u.x[1, 0, 0] - u.x[1, 0, -1]) / Delta) +
                    sq((u.y[1, 0, 1] - u.y[0, 0, 1]) / Delta + (u.x[1, 0, 1] - u.x[1, 0, 0]) / Delta));
    // xz components
    tdef2 += 0.25 * (sq((u.x[0, 1, 0] - u.x[0, 0, 0]) / Delta + (u.z[0, 1, 0] - u.z[-1, 1, 0]) / Delta) +
                    sq((u.x[0, 0, 0] - u.x[0, -1, 0]) / Delta + (u.z[0, 0, 0] - u.z[-1, 0, 0]) / Delta) +
                    sq((u.x[1, 0, 0] - u.x[1, -1, 0]) / Delta + (u.z[1, 0, 0] - u.z[0, 0, 0]) / Delta) +
                    sq((u.x[1, 1, 0] - u.x[1, 0, 0]) / Delta + (u.z[1, 1, 0] - u.z[0, 1, 0]) / Delta));
    // yz components
    tdef2 += 0.25 * (sq((u.z[0, 0, 1] - u.z[0, 0, 0]) / Delta + (u.y[0, 0, 1] - u.y[0, -1, 1]) / Delta) +
                    sq((u.z[0, 0, 0] - u.z[0, 0, -1]) / Delta + (u.y[0, 0, 0] - u.y[0, -1, 0]) / Delta) +
                    sq((u.z[0, 1, 0] - u.z[0, 1, -1]) / Delta + (u.y[0, 1, 0] - u.y[0, 0, 0]) / Delta) +
                    sq((u.z[0, 1, 1] - u.z[0, 1, 0]) / Delta + (u.y[0, 1, 1] - u.y[0, 0, 1]) / Delta));

    // ------------------------------------------------------------
    // STEP 3: Calculate TKE source/sink terms
    // ------------------------------------------------------------
    // Shear production: P = Km * S²
    // Buoyancy production: B = -Kh * N²
    // Combined: (Cm*λ*S²/2 - 3*Cm*λ*N²/2)
    e12p[] = (Cm * lambda * tdef2 / 2 - 3 * Cm * lambda * dbdz / 2);

    // Dissipation: ε = (ce1 + ce2*λ/Δ) * e^(3/2) / λ
    e12p[] -= (ce1 + ce2 * lambda / Delta) * sq(e120[]) / (2 * lambda);

    // Canopy drag sink (wake production of turbulence)
    #if CANOPY
    e12p[] = e12p[] - 4. / 3. * Cd * PAD(y) * U[] * e120[] * CL[];
    #endif
  }

  boundary({Evis});

  // ------------------------------------------------------------
  // STEP 4: Interpolate diffusivities to faces
  // ------------------------------------------------------------
  foreach_face() {
    Km.x[] = (Evis[] + Evis[-1]) / 2;     // Eddy viscosity
    Kh.x[] = Km.x[] / Pr.x[];             // Eddy diffusivity (Pr=1 typically)
    Ke.x[] = Km.x[] * 2.;                 // TKE diffusivity (larger than momentum)
  }

  boundary({Km, Kh, Ke});
}

// ============================================================
// EVENT: tracer_diffusion - Solve TKE diffusion equation
// ============================================================
mgstats mgb;                             // Multigrid solver statistics

event tracer_diffusion(i++)
{
  // Solve TKE evolution: ∂e/∂t + u·∇e = ∇·(Ke∇e) + P + B - ε - Dc
  mgb = diffusion(e120, dt, Ke, r = e12p);
}
