# Part 3: CO2 Transfer to a Leaf - Flow and Scalar Transport

In this tutorial, we mainly want to investigate the influence of wind and leaf orientation on the CO2 diffusive flux to an elliptical leaf. Please do note this is not photosynthetic rate A. We here focus on the maximum CO2 diffusive flux to the leaf surface.

We will couple flow and scalar transport of Basilisk's Navier-Stokes solver with embedded boundary method. The embedded boundary method allows us to handle the leaf geometry as a solid object in the fluid domain.

---

## Overview

In this tutorial, you will learn:

1. **Navier-Stokes Flow**: Solve incompressible flow around complex geometries
2. **Embedded Boundary Method**: Handle solid objects (leaf) in fluid domains
3. **Scalar Transport**: Model advection-diffusion of CO2 in moving air
4. **Surface Flux Calculation**: Compute diffusive flux at embedded boundaries
5. **Parameter Studies**: Investigate effects of wind and leaf orientation

---

## Problem Setup

### General Description

We simulate CO2 diffusive flux from the ambient air to an elliptical leaf surface. The leaf is represented as an embedded solid boundary in a 2D fluid domain. The simulation solves coupled Navier-Stokes and scalar transport equations to capture the interaction between flow and concentration fields.

**Mathematical Framework**:

```
Navier-Stokes:    ∂u/∂t + u·∇u = -∇p + (1/Re)∇²u
                  ∇·u = 0
Scalar transport: ∂s/∂t + u·∇s = ∇·(D∇s)
Boundary:         u = 0 at leaf surface (no-slip)
                  s = s_ls at leaf surface (Dirichlet)
```

### Common Parameters (All Exercises)

| Parameter | Value | Description |
|-----------|-------|-------------|
| Domain size | 120 × 120 units | Computational domain, centered at origin |
| Reynolds number | 50 | Re = U·L/ν (laminar flow regime) |
| Ambient CO2 | 40 mmol/m³ | CO2 concentration in incoming air |
| Leaf surface CO2 | 20 mmol/m³ | CO2 concentration at leaf surface |
| Max grid level | 9 | Maximum adaptive mesh refinement level |
| Simulation time | 0 to 400 | Time units |
| Flux output | t = 300 to 400 | CO2 flux recorded every time unit |

### Key Physics

1. **CO2 Concentration Gradient**: The difference between ambient (40 mmol/m³) and leaf surface (20 mmol/m³) drives diffusive flux into the leaf
2. **Flow Effect**: Wind enhances mass transfer by reducing boundary layer thickness
3. **Geometry Effect**: Leaf orientation affects the interaction with flow field
4. **Embedded Boundary**: The leaf is treated as a solid obstacle with no-slip conditions

---

## Exercise Structure: 2×2 Parameter Study

We conduct four simulations varying two parameters:
- **Wind**: No wind (0 m/s) vs. With wind (1 m/s)
- **Leaf Orientation**: Horizontal vs. Vertical

### Parameter Matrix

| Exercise | Wind Speed | Leaf Orientation | Major Axis | Minor Axis | File Name |
|----------|-----------|------------------|------------|------------|-----------|
| **Ex1**  | 0 m/s (No wind) | Horizontal | r1 = 5 (x-direction) | r2 = 1 | `Ex1_NW_Hleaf.c` |
| **Ex2**  | 1 m/s (With wind) | Horizontal | r1 = 5 (x-direction) | r2 = 1 | `Ex2_W_Hleaf.c` |
| **Ex3**  | 0 m/s (No wind) | Vertical | r1 = 1 | r2 = 5 (y-direction) | `Ex3_NW_Vleaf.c` |
| **Ex4**  | 1 m/s (With wind) | Vertical | r1 = 1 | r2 = 5 (y-direction) | `Ex4_W_Vleaf.c` |

### How to Modify Parameters Between Exercises

The exercises differ only in two key parameters in the source code:

#### 1. Wind Speed (Parameter: `wind_in`)

```c
// Ex1 & Ex3: No wind
double wind_in = 0.;        // Still air condition

// Ex2 & Ex4: With wind
double wind_in = 1.;        // Horizontal wind at 1 m/s
```

**Location in code**: Line ~37 in the "SIMULATION PARAMETERS" section

#### 2. Leaf Orientation (Parameters: `r1`, `r2`)

```c
// Ex1 & Ex2: Horizontal leaf (major axis along x)
double r1 = 5, r2 = 1;      // Ellipse: wide in x, narrow in y

// Ex3 & Ex4: Vertical leaf (major axis along y)
double r1 = 1, r2 = 5;      // Ellipse: narrow in x, wide in y
```

**Location in code**: Line ~53 in the "LEAF GEOMETRY" section

**Note**: The ellipse is defined by the level-set function:
```c
#define ELLIPSE (sq(x/r1) + sq(y/r2) - 1.)
```
Where `sq(x/r1)` means (x/r1)². The leaf surface is where ELLIPSE = 0.

---

## Running the Exercises

### Run Individual Exercises

```bash
# Exercise 1: No wind, Horizontal leaf
make Ex1_NW_Hleaf.tst

# Exercise 2: With wind, Horizontal leaf
make Ex2_W_Hleaf.tst

# Exercise 3: No wind, Vertical leaf
make Ex3_NW_Vleaf.tst

# Exercise 4: With wind, Vertical leaf
make Ex4_W_Vleaf.tst
```

### Run All Exercises

```bash
make Ex1_NW_Hleaf.tst Ex2_W_Hleaf.tst Ex3_NW_Vleaf.tst Ex4_W_Vleaf.tst
```

### Expected Output for Each Exercise

Each exercise produces:

1. **Visualization**: `s.mp4` in the exercise folder
   - Shows CO2 concentration field evolution
   - Color map: blue (20 mmol/m³) to red (40 mmol/m³)
   - Leaf interior is masked

2. **Flux Data**: `diag1` in the exercise folder
   - Two columns: time, CO2_flux
   - Records from t=300 to t=400
   - Used for quantitative comparison

---

## Analyzing and Comparing Results

### Visualize All Results Together

Use the provided Python script to compare CO2 fluxes across all four exercises:

```bash
python plot_diag1.py
```

**The script will**:
- Read `diag1` files from Ex1, Ex2, Ex3, Ex4 folders
- Plot CO2 flux vs. time for all cases
- Save comparison plot as `co2_flux_comparison.png`

### Key Questions to Investigate

1. **Effect of Wind** (Compare Ex1 vs Ex2, Ex3 vs Ex4):
   - How much does wind increase CO2 flux?
   - Does the flux reach steady state faster with wind?
   - What happens to the concentration boundary layer?

2. **Effect of Leaf Orientation** (Compare Ex1 vs Ex3, Ex2 vs Ex4):
   - Does orientation matter without wind?
   - Does orientation matter with wind?
   - Which orientation is better for CO2 uptake?

3. **Combined Effects** (Compare all four):
   - Which configuration gives the highest CO2 flux?
   - Is there an interaction between wind and orientation?
   - What is the optimal condition for maximum gas exchange?

### Expected Trends

Based on mass transfer theory, we expect:

- **Wind increases flux**: Thinner boundary layer → steeper concentration gradient → higher flux
- **Orientation matters with wind**: Vertical leaf has larger cross-sectional area perpendicular to flow
- **Orientation doesn't matter without wind**: Pure diffusion is isotropic
- **Highest flux**: Ex4 (With wind + Vertical leaf) should give maximum flux

**Your task**: Verify these expectations with simulation results!

---

## Bonus Challenges

Ready for more? We've prepared **5 bonus challenges** with complete solutions in the `bonus/` folder!

### Quick Start

```bash
cd bonus
make help           # See all available challenges
make Bonus1_VaryRe.tst  # Run a specific challenge
```

### Challenge Overview

| # | Challenge | File | Key Variable |
|---|-----------|------|--------------|
| 1 | **Vary Reynolds Number** | `Bonus1_VaryRe.c` | Re = 10, 50, 100, 500 |
| 2 | **Circular Leaf Shape** | `Bonus2_CircularLeaf.c` | Circle vs. Ellipse |
| 3 | **Vary Concentration** | `Bonus3_VaryConcentration.c` | Different ΔC values |
| 4 | **Inclined Leaf** | `Bonus4_InclinedLeaf.c` | Rotation angle θ |
| 5 | **Multiple Leaves** | `Bonus5_MultipleLeaves.c` | Leaf spacing |

### Detailed Descriptions

1. **Vary Reynolds number**: How does flow regime affect flux? Test Re = 10, 50, 100, 500 to see the transition from viscous to inertial flow.

2. **Different leaf shapes**: Is the ellipse more efficient than a circle? Compare CO2 uptake for leaves of the same surface area.

3. **Vary concentration gradient**: Test Fick's law! Is flux linear with ΔC? Try different ambient and leaf surface concentrations.

4. **Inclined leaf**: Find the optimal angle! Rotate the ellipse by 0°, 30°, 45°, 60°, 90° to find which gives maximum flux.

5. **Multiple leaves**: Study interaction effects! How does a downstream leaf perform in the wake of an upstream leaf?

### More Information

See `bonus/README.md` for:
- Detailed instructions for each challenge
- Expected results and analysis tasks
- Tips for parameter sweeps and advanced analysis
- Plotting and visualization suggestions

---

## Understanding the Physics

### 1. Embedded Boundary Method

The leaf is represented using a level-set function:

```c
#define ELLIPSE (sq(x/r1) + sq(y/r2) - 1.)
```

- **ELLIPSE < 0**: Inside the leaf (solid)
- **ELLIPSE = 0**: Leaf surface
- **ELLIPSE > 0**: Outside the leaf (fluid)

Basilisk computes volume fractions (`cs`) and face fractions (`fs`) from this level-set to handle the boundary accurately.

### 2. Boundary Conditions at Leaf Surface

```c
u.t[embed] = dirichlet(0.);   // No tangential velocity (no-slip)
u.n[embed] = dirichlet(0.);   // No normal velocity (no penetration)
s[embed] = dirichlet(s_ls);   // Fixed CO2 concentration = 20 mmol/m³
```

These conditions ensure:
- The leaf is a solid obstacle to flow
- CO2 concentration at leaf surface is maintained at 20 mmol/m³
- Diffusive flux is driven by concentration gradient at the surface

### 3. CO2 Flux Calculation

The total CO2 flux into the leaf is computed using:

```c
event diag_flux (t = 300; t += 1) {
  double flx = 0;
  foreach(reduction(+:flx)) {
    double val = 0, e = embed_flux(point, s, mu, &val);
    if (val)
      flx += (val - e*s[])*sq(Delta);  // Sum diffusive flux
  }
  fprintf(fp, "%g %g\n", t, flx);       // Write to diag1 file
}
```

This calculates: **Flux = -D ∇s · n** integrated over the entire leaf surface.

Where:
- **D**: Effective diffusivity (scaled by viscosity/Re)
- **∇s**: CO2 concentration gradient at the surface
- **n**: Surface normal vector pointing into the leaf

### 4. Adaptive Mesh Refinement

The grid is dynamically refined based on solution features:

```c
adapt_wavelet({cs, s, u}, (double[]){5e-1, be, ue, ue}, maxlevel, 5);
```

Refinement occurs where:
- **cs** has high gradients (near leaf surface)
- **s** has high gradients (concentration boundary layer)
- **u** has high gradients (velocity boundary layer)

This focuses computational effort where it's needed most.

---

## Boundary Layer Theory

### Concentration Boundary Layer

The concentration boundary layer thickness (δ_c) scales as:

```
δ_c ~ L / √(Re·Sc)
```

Where:
- **L**: Leaf characteristic length (~5)
- **Re**: Reynolds number (50)
- **Sc**: Schmidt number (ν/D ~ 1 in this simulation)

**Key insight**: Higher wind → higher Re → thinner boundary layer → steeper gradient → higher flux

### Flow Patterns

**Without wind** (Ex1, Ex3):
- Symmetric concentration field around leaf
- Pure diffusion dominates
- Slow CO2 transport

**With wind** (Ex2, Ex4):
- Asymmetric concentration field (high upstream, low downstream)
- Advection dominates near leaf
- Fast CO2 transport
- Wake formation downstream

---
