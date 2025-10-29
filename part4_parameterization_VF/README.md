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

   - `dump-080`: Simulation snapshot at t=80s (for restart)
   - `W12/slice_80`: 2D slice data at end time
3. **Visualization**:

   You can plot vertical profiles and 2D slices of the fields using the Python script:

   ```bash
   # Activate conda environment
   conda activate workshop

   # Run plotting script
   python plot_slice_profiles.py
   ```

   This will generate:
   - Vertical profiles of `cw` and `ux`
   - 2D slices showing spatial distribution (with roof masking)
   - Output files: `profile_cw.png`, `profile_Ux.png`, `2D_slice_cw.png`, `2D_slice_ux.png`

4. **Directory structure**:

   ```
   part4_parameterization_VF/
   ├── Green2D.c              # Main simulation
   ├── physics.h              # Physical parameterizations
   ├── Canopy.h               # Vegetation model
   ├── SGS_TKE.h              # Turbulence model
   ├── Makefile               # Build instructions
   ├── plot_slice_profiles.py # Visualization script
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

## Exercise: Wind Effect on Canopy Flow

### **Exercise 1: Test Different Wind Speeds**

**Objective:** Understand how background wind speed affects flow patterns around the canopy.

**Step 1: Create experiment directories**

```bash
# Create main experiment folder
mkdir wind_exp
cd wind_exp

# Create two sub-folders for different wind conditions
mkdir W1 W2

# Copy source files to each folder
cp ../*.c ../*.h W1/
cp ../*.c ../*.h W2/
cp ../Makefile W1/
cp ../Makefile W2/
```

**Step 2: Modify wind conditions**

Edit the `physics.h` file in each folder to change the wind speed (U0 parameter, line 27):

**For W1 folder (Weak wind):**
```bash
cd W1
# Edit physics.h, change line 27:
#define U0 0.2    // Weak wind: 0.2 m/s
```

**For W2 folder (Strong wind):**
```bash
cd W2
# Edit physics.h, change line 27:
#define U0 0.6    // Strong wind: 0.6 m/s
```

**Step 3: Run simulations**

```bash
# Run W1 simulation
cd W1
make Green2D.tst
cd ..

# Run W2 simulation
cd W2
make Green2D.tst
cd ..
```

**Step 4: Visualize and compare**

```bash
conda activate workshop

# Plot W1 results
cd W1
python ../plot_slice_profiles.py
cd ..

# Plot W2 results
cd W2
python ../plot_slice_profiles.py
cd ..
```

**Questions to consider:**
1. How does wind speed affect the velocity field around the canopy?
2. Where do you see wake regions forming behind vegetation?
3. How does stronger wind change the water vapor distribution?
4. Compare the 2D slices - what differences do you observe?

**Expected learning:** Background wind creates wakes behind obstacles. Stronger wind leads to more pronounced wake regions and affects scalar transport.

---

## Additional Exercises to Explore

Feel free to explore other parameters on your own:

1. **Canopy Geometry**: Modify canopy positions (CUBE1_X, CUBE2_X, CUBE3_X) or dimensions (CUBE_HEIGHT, CUBE_WIDTH) in `Canopy.h`
2. **Roof Shape**: Change NUM_WAVES or ROOF_AMPLITUDE in `Green2D.c` to see how roof geometry affects circulation
3. **Heating**: Increase leaf temperature (TV[] initialization) or surface buoyancy (BSURF) to observe thermal convection
4. **Boundary Conditions**: Switch from periodic to inflow/outflow boundaries in `physics.h`

Experiment with these parameters and observe how they change the greenhouse microclimate!

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
