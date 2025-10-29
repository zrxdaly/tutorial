# Part 4: Greenhouse Flow with Canopy Parameterization

In this tutorial, we explore indoor climate simulation in a greenhouse with vegetation canopy and a wave-shaped roof. This simulation demonstrates how to couple multiple physical processes: turbulent flow, canopy energy balance, scalar transport (water vapor), and geometric constraints (roof boundary).

We build upon concepts from previous tutorials and introduce **sub-grid scale turbulence modeling** and **vegetation parameterization** for Large Eddy Simulation (LES).

---

## Overview of the Code Structure

This simulation is organized into four main files, each handling specific physical processes:

### 1. **Green2D.c** - Main Simulation Driver

This is the "control center" of your simulation. Think of it as the conductor of an orchestra.

**What it does:**

- Sets up the computational domain (100m × 100m greenhouse)
- Defines the wave-shaped roof boundary that blocks flow
- Controls when different physical processes are computed
- Generates visualization outputs (videos of buoyancy, velocity, water vapor)
- Saves simulation snapshots for later analysis

**Key features:**

- **Roof geometry**: A wavy ceiling (like a sine wave) from y=70m to y=100m that forces velocity to zero
- **Adaptive mesh**: Grid automatically refines near important features (canopy, roof)
- **Event-driven**: Different physics happen at different times (e.g., roof damping every timestep, videos every second)

**For beginners:** Start here to understand the simulation setup and what outputs you'll get.

---

### 2. **physics.h** - Physical Parameterizations

This file handles the "rules of physics" for the indoor climate.

**What it does:**

- Defines physical constants (gravity, air properties, initial wind)
- Sets boundary conditions (periodic sides, no-slip top/bottom)
- Computes buoyancy force (warm air rises, cool air sinks)
- Applies canopy drag (vegetation slows down the wind)
- Handles scalar diffusion with source terms (heat and water vapor from plants)

**Key physics:**

- **Buoyancy**: Converts temperature differences to vertical acceleration
- **Canopy drag**: Trees act like obstacles that slow the wind (drag = Cd × PAD × |u| × u)
- **Heat and moisture sources**: Plants release heat and water vapor into the air

**For beginners:** This connects the math (equations) to the simulation. If you want to change wind speed or add heating, modify this file.

---

### 3. **Canopy.h** - Vegetation Energy Balance

This file treats vegetation as active participants in the flow, not just obstacles.

**What it does:**

- Defines three cube-shaped canopy elements (like trees or crop rows)
- Computes leaf temperature by balancing multiple energy fluxes
- Calculates heat exchange between leaves and air
- Computes transpiration (water vapor released by plants)

**Key physics (3-step energy balance):**

1. **Geometry**: Uses the "fractions method" to define where vegetation is located
2. **Radiation**: Longwave radiation exchange between leaves, sky, and ground
3. **Convection & Transpiration**:
   - Heat flux: H = (T_leaf - T_air) / resistance
   - Water flux: E = (humidity_saturation - humidity_air) / (resistance + stomatal_resistance)

**For beginners:** Plants are not passive! They heat up from sunlight and cool down by releasing water vapor (like sweating). This file models that process.

---

### 4. **SGS_TKE.h** - Sub-Grid Scale Turbulence

This file handles turbulence that is too small for the grid to resolve directly.

**What it does:**

- Computes turbulent kinetic energy (TKE = energy of swirling eddies)
- Calculates "eddy viscosity" (how much turbulence mixes momentum)
- Determines mixing length (how far turbulent eddies can reach)
- Accounts for turbulence production (from shear and buoyancy) and dissipation

**Key concepts:**

- **Mixing length**: In unstable air (warm below), eddies can grow large. In stable air (warm above), they're suppressed
- **TKE equation**: ∂e/∂t = Shear Production + Buoyancy Production - Dissipation - Canopy Drag
- **Stability-dependent**: The code adjusts turbulence based on whether air is stable or unstable

**For beginners:** Turbulence is complex! This file uses a simplified model: instead of resolving every tiny eddy, we track the average energy of eddies and use that to estimate mixing.

---

## How the Files Work Together

```
Green2D.c (main)
    ↓
    Includes physics.h
        ↓
        Includes Canopy.h  →  Defines vegetation geometry & energy balance
        ↓
        Includes SGS_TKE.h  →  Computes turbulent mixing
    ↓
    Each timestep:
        1. SGS_TKE.h: Update eddy viscosity & TKE
        2. Canopy.h: Compute leaf energy balance
        3. physics.h: Apply buoyancy, drag, and diffusion
        4. Green2D.c: Apply roof boundary condition
        5. Adapt grid based on solution features
```

---

## Running the Simulation

### Quick Start

```bash
# Compile and run with 4 MPI processes
make Green2D.tst

# Clean up output files
make cleanout

# Clean everything (compiled files + outputs)
make cleanall
```

### Expected Outputs

After running, you'll get:

1. **Videos** (every 1 second of simulation time):

   - `b.mp4`: Buoyancy field (shows temperature stratification)
   - `ux.mp4`: Horizontal velocity (shows wind patterns)
   - `cw.mp4`: Water vapor concentration
2. **Data files**:

   - `dump-00005`: Simulation snapshot at t=5s (for restart)
   - `W12/slice_5`: 2D slice data at end time
3. **Directory structure**:

   ```
   part4_parameterization_VF/
   ├── Green2D.c              # Main simulation
   ├── physics.h              # Physical parameterizations
   ├── Canopy.h               # Vegetation model
   ├── SGS_TKE.h              # Turbulence model
   ├── Makefile               # Build instructions
   ├── b.mp4, ux.mp4, cw.mp4  # Visualization outputs
   └── W12/                   # Output directory
   ```

---

## Understanding the Simulation Setup

### Domain Configuration

| Parameter       | Value                  | Description                               |
| --------------- | ---------------------- | ----------------------------------------- |
| Domain size     | 100m × 100m           | Horizontal × Vertical                    |
| Grid resolution | 64 base + 7 levels AMR | ~128×128 effective near features         |
| Time duration   | 5 seconds              | Short for demonstration                   |
| Physics         | Buoyancy-driven flow   | Temperature differences drive circulation |

### Canopy Configuration

- **Three cube elements**: Located at x = 20m, 50m, 80m
- **Size**: 5m wide × 30m tall (like trees or tall crops)
- **Properties**: Release heat and water vapor based on energy balance
- **Effect**: Slow down wind, add turbulence, modify temperature and humidity

### Roof Configuration

- **Shape**: Wavy roof with 2 complete waves across domain
- **Height**: Varies from 70m (lowest) to 100m (highest)
- **Effect**: Forces velocity to zero above the wave (solid boundary)
- **Purpose**: Simulate enclosed greenhouse environment

---

## Interesting Exercises for Students

### 🌱 **Exercise 1: Modify Canopy Geometry**

**Objective:** Understand how vegetation distribution affects flow patterns.

**Tasks:**

1. **Change canopy positions** (Canopy.h, lines 26-28):

   ```c
   #define CUBE1_X 30.   // Move first cube
   #define CUBE2_X 50.   // Keep second cube
   #define CUBE3_X 70.   // Move third cube closer
   ```
2. **Change canopy size** (Canopy.h, lines 21-22):

   ```c
   #define CUBE_HEIGHT 20.   // Shorter vegetation
   #define CUBE_WIDTH 10.    // Wider vegetation
   ```
3. **Run simulation and compare**:

   - How does velocity change between vegetation elements?
   - Where do you see the strongest turbulence?
   - How does spacing affect the flow?

**Expected learning:** Vegetation acts as obstacles. Closer spacing creates more drag and turbulence.

---

### 🌊 **Exercise 2: Modify Roof Shape**

**Objective:** Explore how roof geometry affects flow circulation.

**Tasks:**

1. **Change wave number** (Green2D.c, line 28):

   ```c
   #define NUM_WAVES 1    // Single wave (gentle)
   #define NUM_WAVES 4    // Four waves (complex)
   ```
2. **Change wave amplitude** (Green2D.c, line 32):

   ```c
   #define ROOF_AMPLITUDE 5.    // Flatter roof
   #define ROOF_AMPLITUDE 20.   // Steeper roof
   ```
3. **Visualize and analyze**:

   - How does flow pattern change with roof shape?
   - Does a wavy roof create recirculation zones?
   - Where does air get "trapped"?

**Expected learning:** Geometry strongly influences flow. Complex shapes create complex flow patterns.

---

### 🌡️ **Exercise 3: Heat the Canopy**

**Objective:** Observe buoyancy-driven convection from heated vegetation.

**Tasks:**

1. **Increase leaf temperature** (Green2D.c, line 109):

   ```c
   TV[] = 305.15;   // Warm leaves (32°C instead of 22°C)
   ```
2. **Increase ground temperature** (physics.h, line 30):

   ```c
   #define BSURF (gCONST / TREF * 30.0)  // Hotter surface
   ```
3. **Run and observe**:

   - Watch `b.mp4`: Do you see rising plumes above vegetation?
   - Check `ux.mp4`: Does buoyancy create horizontal flow?
   - How far do thermal plumes reach before hitting the roof?

**Expected learning:** Heating creates convection. Warm air rises, creating circulation patterns.

---

### 💨 **Exercise 4: Add Background Wind**

**Objective:** See how mean wind interacts with canopy and roof.

**Tasks:**

1. **Add initial wind** (physics.h, line 28):

   ```c
   #define WIND(s) (2.0)   // Constant 2 m/s wind
   ```
2. **Change boundary conditions** (physics.h, uncomment lines 77-80):

   ```c
   u.n[left] = dirichlet(2.0);    // Inflow at left
   p[left] = neumann(0.);
   u.n[right] = neumann(0.);      // Outflow at right
   p[right] = dirichlet(0.);
   ```

   And comment out `periodic(left);` on line 74
3. **Analyze**:

   - How does wind interact with vegetation?
   - Do you see wake regions behind each canopy element?
   - How does the wavy roof deflect the wind?

**Expected learning:** Wind + obstacles = wakes. Flow separates behind vegetation creating complex patterns.

---

## Physical Insights

### What This Simulation Teaches

1. **Coupled Physics**: Real systems involve many interacting processes (flow, heat, moisture, turbulence)
2. **Scale Separation**:

   - Large scales (mean flow, canopy spacing): Resolved explicitly
   - Small scales (turbulent eddies): Modeled via SGS_TKE
3. **Parameterization Philosophy**:

   - We can't simulate every leaf, so we model vegetation as a porous medium
   - We can't resolve every eddy, so we model turbulence statistically
4. **Geometry Matters**:

   - Canopy distribution affects mixing
   - Roof shape creates recirculation
   - Both influence greenhouse microclimate
5. **Engineering Applications**:

   - Greenhouse design for optimal crop growth
   - Urban canopy flows (trees in cities)
   - Indoor air quality and ventilation

---

## Further Reading

### Key Papers

- **Canopy flows:** Patton et al. (2016), Boekee et al. (2023)
- **SGS turbulence:** Deardorff (1980), Dai et al. (2014)
- **Embedded boundaries:** Popinet (2009) - Basilisk documentation

### Concepts to Explore

- Large Eddy Simulation (LES)
- Plant Area Density (PAD)
- Stomatal conductance and transpiration

---

## Summary

This simulation demonstrates a complete LES framework for greenhouse indoor climate:

- Resolved flow field with adaptive mesh refinement
- Sub-grid scale turbulence closure
- Canopy energy balance and source terms
- Complex geometry (wavy roof)
- Coupled heat and moisture transport

**Key takeaway:** Indoor climate flows are complex, but by breaking the problem into modular components (main driver, physics, canopy, turbulence), we can build comprehensive simulations that capture the essential processes.
