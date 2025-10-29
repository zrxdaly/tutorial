# Part 2: Diffusion Modeling in Basilisk

This tutorial introduces diffusion processes using both manual finite difference methods and Basilisk's built-in solvers, demonstrated through CO2 dispersion modeling.

---

## Overview

In this tutorial, you will learn:

1. **Finite Difference Method**: Manually implement the 2D diffusion equation using central differences
2. **Basilisk Diffusion Solver**: Use Basilisk's implicit solver for stable and efficient diffusion
3. **Source Terms**: Model injection sources in diffusion problems
4. **Dynamic Visualization**: Create animations of time-evolving concentration fields

---

# Tutorial Structure

### Exercise 1: Manual Finite Difference Implementation (`Ex1_FDCO2.c`)

**Objective**: Implement CO2 diffusion using explicit finite difference method.

**Key Concepts**:

- Explicit time integration for diffusion equation: ∂C/∂t = D∇²C
- Central difference approximation for second derivatives
- Stability condition: Δt ≤ Δx²/(4D)
- Source term implementation (CO2 injection)

**Physical Setup**:

- Domain: 5m × 5m
- Initial condition: 400 ppm CO2 (atmospheric background)
- Injection source: Pipe at center (radius 0.1m) maintaining 1500 ppm
- Diffusivity coefficient: D = 1.60e-5 m²/s

**Mathematical Implementation**:

```c
// Second derivatives using central differences
d²C/dx² ≈ (C[i+1,j] - 2C[i,j] + C[i-1,j]) / Δx²
d²C/dy² ≈ (C[i,j+1] - 2C[i,j] + C[i,j-1]) / Δy²

// Time update
C_new = C_old + Δt × D × (d²C/dx² + d²C/dy²)
```

**Run the example**:

```bash
make Ex1_FDCO2.tst
```

**Question**:
Is the simulation time long enough to let the CO2 diffuse from the injection point to the edges of the room? is that realistic?

**Output**: `CO2_evo.mp4` showing CO2 spreading from injection point

---

### Exercise 2: Basilisk Diffusion Solver (`Ex2_basi.c`)

**Objective**: Solve the same diffusion problem using Basilisk's implicit solver.

**Key Concepts**:

- Implicit solver: More stable, allows larger timesteps
- `diffusion.h` solver usage
- Face vector for diffusion coefficients

**Key Advantages of Implicit Solver**:

- **mostly stable**: No strict timestep limitation
- **Efficient**: Built-in sparse matrix solvers

**Solver Usage**:

```c
const face vector kappa[] = {D, D};  // Diffusion coefficient in x,y
diffusion(C, dt, kappa);             // Solve: dC/dt = D∇²C
```

**Run the example**:

```bash
make Ex2_basi.tst
```

**Compare**:

- Same physical setup as Exercise 1
- Different numerical method
- Observe stability and accuracy differences
- an example of using external solvers

**Output**: `CO2_field.mp4` showing CO2 concentration evolution

---

### Exercise 3: Multiple Injection Sources (`Ex3_inject3.c`)

**Objective**: Extend the model to include three CO2 injection points.

**Key Concepts**:

- Multiple source terms
- Source interaction and superposition
- Spatial distribution of sources
- Data output for quantitative analysis

**Source Configuration**:

- Nozzle 1: Center (0, 0)
- Nozzle 2: Left (-1.25, 0)
- Nozzle 3: Right (1.25, 0)
- All maintain 1500 ppm within 0.1m radius

**New Feature - Data Export**:

```c
event printdata (t = 0; t <= 10.; t += 0.3) {
  // Output concentration along centerline for analysis
  for (double x = -L0/2; x < L0/2; x += L0/N)
    fprintf(fp, "%g %g %g\n", t, x, interpolate(C, x, 0));
}
```

**Run the example**:

```bash
make Ex3_inject3.tst
```

**Output**:

- `CO2_field.mp4`: Visualization of three plumes merging
- `y0.dat`: Centerline concentration data for plotting


### Bonus Questions

1. How does the spacing between nozzles affect:

   - Time to reach target concentration?
   - Uniformity of distribution?
   - Peak concentration at the center?
2. Can you achieve >80% coverage (800+ ppm) within 3 seconds? What configuration works best?
3. If you could only use 2 nozzles instead of 3, where would you place them?

### Suggested Configurations to Try

1. **Linear arrangement**: All three along x-axis at different spacings
2. **Triangular pattern**: Equilateral triangle formation
3. **L-shape**: Corner coverage strategy
4. **Your creative design**: What works best?

---
