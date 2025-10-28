# Bonus Challenges - Advanced Leaf Gas Exchange Simulations

This folder contains solutions to five bonus challenges that extend the basic leaf CO2 flux simulations with more advanced scenarios and parameter studies.

---

## Overview of Bonus Challenges

| Challenge | File | Main Variable | Key Question |
|-----------|------|---------------|--------------|
| **Bonus 1** | `Bonus1_VaryRe.c` | Reynolds number (Re) | How does flow regime affect flux? |
| **Bonus 2** | `Bonus2_CircularLeaf.c` | Leaf shape | Is ellipse better than circle? |
| **Bonus 3** | `Bonus3_VaryConcentration.c` | Concentration gradient | Is flux linear with ΔC? |
| **Bonus 4** | `Bonus4_InclinedLeaf.c` | Leaf inclination angle | What angle gives maximum flux? |
| **Bonus 5** | `Bonus5_MultipleLeaves.c` | Leaf spacing | Do leaves interact? |

---

## Bonus Challenge 1: Vary Reynolds Number

**File**: `Bonus1_VaryRe.c`

**Objective**: Investigate how flow regime (Re) affects CO2 flux and the flux-orientation relationship.

### What to Modify

In the "BONUS CHALLENGE" section (around line 30), uncomment different Re values:

```c
double Re = 10;      // Very viscous, thick boundary layer
// double Re = 50;      // Moderate (default case)
// double Re = 100;     // Transitional
// double Re = 500;     // More inertial, thin boundary layer
```

### Run the Challenge

```bash
cd bonus
make Bonus1_VaryRe.tst
```

### Expected Results

| Re | Boundary Layer | Expected Flux | Flow Regime |
|----|----------------|---------------|-------------|
| 10 | Very thick | Lower | Viscous dominated |
| 50 | Moderate | Baseline | Transitional |
| 100 | Thin | Higher | More inertial |
| 500 | Very thin | Highest | Inertial dominated |

**Theoretical prediction**: Flux scales as √Re due to boundary layer thinning.

### Analysis Tasks

1. Run simulations for Re = 10, 50, 100, 500
2. Plot flux vs. Re on log-log scale
3. Check if flux ∝ Re^α (find the exponent α)
4. Compare horizontal vs. vertical leaf at each Re
5. Does orientation matter more at high or low Re?

---

## Bonus Challenge 2: Circular Leaf

**File**: `Bonus2_CircularLeaf.c`

**Objective**: Compare CO2 flux between circular and elliptical leaves of the same surface area.

### What's Different

The ellipse is replaced with a circle:

```c
// Original ellipse: area = π * 5 * 1 ≈ 15.7
// Circle: r = √5 ≈ 2.236 for same area

double r_circle = 2.236;
#define CIRCLE (sq(x) + sq(y) - sq(r_circle))
```

### Run the Challenge

```bash
cd bonus
make Bonus2_CircularLeaf.tst
```

### Expected Results

- **Without wind**: Circle and ellipse should have similar flux (isotropic diffusion)
- **With wind**: Ellipse (horizontal) should have higher flux (larger projected area)

### Analysis Tasks

1. Run circular leaf with wind_in = 0 and wind_in = 1
2. Compare with Ex1 (horizontal ellipse, no wind) and Ex2 (with wind)
3. Calculate flux per unit area for each shape
4. Which shape is more "efficient"?
5. Does the answer depend on wind speed?

---

## Bonus Challenge 3: Vary Concentration Gradient

**File**: `Bonus3_VaryConcentration.c`

**Objective**: Test whether CO2 flux follows Fick's law (Flux ∝ ΔC).

### What to Modify

In the "BONUS CHALLENGE" section (around line 26), uncomment different concentration cases:

```c
// Case 1: Default gradient (ΔC = 20)
double s_in = 40.0;
double s_ls = 20.0;

// Case 2: Large gradient (ΔC = 50)
// double s_in = 60.0;
// double s_ls = 10.0;

// Case 3: Small gradient (ΔC = 5)
// double s_in = 30.0;
// double s_ls = 25.0;

// Case 4: Very large gradient (ΔC = 80)
// double s_in = 100.0;
// double s_ls = 20.0;
```

### Run the Challenge

```bash
cd bonus
make Bonus3_VaryConcentration.tst
```

### Expected Results

According to Fick's law of diffusion:
```
Flux = -D * ∇C
```

For a given geometry, flux should be **linear** with concentration difference (ΔC = s_in - s_ls).

### Analysis Tasks

1. Run simulations for all four cases
2. Plot flux vs. ΔC
3. Fit a linear regression: Flux = k * ΔC
4. Check R² value (should be close to 1)
5. Does linearity hold for very large gradients?

**Note**: The output file `diag1` includes three columns: time, flux, and ΔC for easy plotting.

---

## Bonus Challenge 4: Inclined Leaf

**File**: `Bonus4_InclinedLeaf.c`

**Objective**: Find the optimal leaf inclination angle for maximum CO2 flux.

### What to Modify

In the "BONUS CHALLENGE" section (around line 42), uncomment different angles:

```c
double theta = 0.0;         // 0° - Horizontal
// double theta = M_PI/6;      // 30° inclination
// double theta = M_PI/4;      // 45° inclination
// double theta = M_PI/3;      // 60° inclination
// double theta = M_PI/2;      // 90° - Vertical
```

### Run the Challenge

```bash
cd bonus
make Bonus4_InclinedLeaf.tst
```

### Expected Results

- **0° (horizontal)**: Baseline from Ex2
- **45°**: Possibly optimal (compromise between projected area and wake)
- **90° (vertical)**: Same as Ex4

**Hypothesis**: There should be an optimal angle between 0° and 90° that maximizes flux.

### Analysis Tasks

1. Run simulations for θ = 0°, 30°, 45°, 60°, 90°
2. Plot flux vs. angle
3. Find the angle that gives maximum flux
4. Does the optimal angle depend on Re?
5. Explain the result using boundary layer theory

**Advanced**: Run a sweep from 0° to 90° in 10° increments to find the exact optimum.

---

## Bonus Challenge 5: Multiple Leaves

**File**: `Bonus5_MultipleLeaves.c`

**Objective**: Study interaction effects between multiple leaves in the same flow.

### What's Different

Two leaves are defined:

```c
// Leaf 1 (upstream): centered at x = -10
double x1 = -10.0;
double y1 = 0.0;

// Leaf 2 (downstream): centered at x = +10
double x2 = 10.0;
double y2 = 0.0;
```

### What to Modify

Adjust leaf positions to study different configurations:

```c
// Try different spacings:
// double x2 = 8.0;   // Close spacing (spacing = 18)
// double x2 = 15.0;  // Wide spacing (spacing = 25)
// double y2 = 5.0;   // Offset vertically
```

### Run the Challenge

```bash
cd bonus
make Bonus5_MultipleLeaves.tst
```

### Expected Results

- **Upstream leaf (Leaf 1)**: Should have similar flux to single leaf
- **Downstream leaf (Leaf 2)**: Should have **reduced** flux due to:
  1. CO2 depletion by upstream leaf
  2. Wake effects from upstream leaf
  3. Flow disturbance

**Key question**: Is `flux_total < 2 * flux_single_leaf`?

### Analysis Tasks

1. Run single leaf simulation (Ex2) for reference
2. Run double leaf with spacing = 20
3. Compare: `flux_leaf1 + flux_leaf2` vs. `2 * flux_Ex2`
4. Calculate efficiency: `(flux_leaf1 + flux_leaf2) / (2 * flux_Ex2)`
5. Vary spacing: does wider spacing reduce interaction?
6. Try vertical offset (y2 ≠ 0): does it help?

**Note**: The output file `diag1` has four columns: time, total_flux, flux_leaf1, flux_leaf2.

---

## General Instructions

### Compiling and Running

Each bonus challenge can be compiled and run individually:

```bash
cd bonus

# Compile a specific challenge
make Bonus1_VaryRe.tst

# Or run all bonuses
make Bonus1_VaryRe.tst Bonus2_CircularLeaf.tst Bonus3_VaryConcentration.tst \
     Bonus4_InclinedLeaf.tst Bonus5_MultipleLeaves.tst
```

### Output Files

Each simulation produces:

1. **s.mp4**: Video visualization of CO2 concentration field
2. **diag1**: Time series data file with flux measurements

### Plotting Results

You can adapt the `plot_diag1.py` script from the parent directory:

```bash
# Copy the plotting script
cp ../plot_diag1.py .

# Modify it to read bonus folder outputs
# Then run:
python plot_diag1.py
```
