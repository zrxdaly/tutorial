# Part 1: Basic Grid Operations in Basilisk

This tutorial introduces fundamental concepts of computational fluid dynamics using Basilisk: grid generation, field allocation, and adaptive mesh refinement.

---

## Overview

In this tutorial, you will learn:
1. **Grid Generation**: Create 2D Cartesian grids
2. **Field Allocation**: Define and initialize scalar fields on the grid
3. **Adaptive Mesh Refinement (AMR)**: Refine grids in regions with high gradients
4. **Field Visualization**: Visualize temperature distributions and grid structures

---

## Tutorial Structure

### Exercise 1: Uniform Grid Generation (`Ex1_grid.c`)

**Objective**: Create a simple 2D uniform grid and visualize it.

**Key Concepts**:
- Domain definition using `size()` and `origin()`
- Uniform grid initialization with `init_grid()`
- Basic visualization with `cells()` and `box()`

**Run the example**:
```bash
make Ex1_grid.tst
```

---

### Exercise 2: Temperature Field Allocation (`Ex2_var.c`)

**Objective**: Allocate a scalar field (temperature) on the grid and visualize it.

**Key Concepts**:
- Scalar field declaration with `scalar temp[]`
- Field initialization using `foreach()` loops
- Spatial functions: exponential decay from a center point
- Field statistics with `statsf()`

**Initial Setup**: Start with 4×4 grid resolution
```c
init_grid(4);  // 4×4 = 16 cells
```

**Temperature Function**:
```
temp(x,y) = 20.0 + 10.0 × exp(-60.0 × r²)
where r = √[(x - x₀)² + (y - y₀)²]
```
This creates a hot spot (30°C) at the center, decaying to 20°C at the edges.

**Run the example**:
```bash
make Ex2_var.tst
```

** Try This**:
- Change `init_grid(4)` to `init_grid(16)` for 16×16 resolution
- Observe how the temperature field becomes better represented
- Compare the min/max/mean statistics between resolutions

**Output**:
- `Temperature_field.png`: Visualization with color map and grid overlay
- Terminal output: Temperature statistics

**Bonus**: Use the Python script to plot temperature profiles:
```bash
conda activate workshop
python temp_lines.py
```

---

### Exercise 3: Adaptive Grid Refinement (`Ex3_Agrid.c`)

**Objective**: Use adaptive mesh refinement to increase resolution where temperature gradients are highest.

**Key Concepts**:
- Adaptive mesh refinement with `refine()`
- Region-based refinement using geometric conditions
- Multi-level grids with different resolutions

**Key Observation**:
By examining the temperature field from Exercise 2, we see that the steepest gradients occur within a circle of radius ≈ 0.3 around the center. This is where we need higher resolution!

**Refinement Strategy**:
```c
int max_level = 5;
double circle_radius = 0.3;
refine(sq(x - circle_x) + sq(y - circle_y) < sq(circle_radius) && level < max_level);
```

This refines cells inside the circle up to level 5, while keeping coarser cells outside.

**Run the example**:
```bash
make Ex3_Agrid.tst
```

**Compare**:
- Uniform 16×16 grid: 256 cells
- Adaptive grid (level 5 inside, level 4 outside): Fewer cells, better accuracy where needed!

**Output**: `Temperature_field.png` showing refined grid near the hot spot

---

## Exercise 4: Your Turn! - Two Hot Spots Challenge (`Ex4_two_spots.c`)

**Objective**: create TWO temperature sources

### The Challenge

Create a temperature field with **two circular hot spots**:
- Hot spot 1: Center at (-0.2, 0.0), temperature peak 35°C, decay constant 60 
- Hot spot 2: Center at (0.2, 0.0), temperature peak 28°C, decay constant 30
- Background temperature: 18°C


### Getting Started

A template file `Ex4_two_spots.c` is copyed from `Ex3_Agrid.c` with some instructions to complete.

### Steps to Complete

1. **Define parameters** fill in parameters of circle 2
2. **Implement adaptive refinement** modify the refine function including circle 2 
3. **Initialize the temperature field** add function of circle 2

### Run Your Solution

```bash
make Ex4_two_spots.tst
```

### Bonus Questions 

**Experiment**:
   - Move the hot spots closer together
   - Change decay constants to make wider or narrower hot spots
   - Try different refinement strategies

1. How many cells does your adaptive grid use compared to a uniform grid of the same finest resolution?
2. What happens to the temperature at (0, 0) as you move the hot spots closer?
3. a refinement strategy that uses even fewer cells?
---

## Tips and Tricks

### Debugging
- Use `fprintf(stderr, ...)` to print debug information
- Check field statistics to verify your initialization
- Visualize early and often!

### Grid Refinement
- Start with lower max_level (4-5) for faster testing
- Too much refinement can slow down simulations
- Refine where you need it, not everywhere!

### Temperature Functions
- Exponential decay: Sharp transitions, good for hot spots
- Gaussian-like: `exp(-k × r²)` where larger k = faster decay
- Linear: Simple but creates sharp boundaries

### Visualization
- Use `cells()` to see the grid structure
- Use `squares()` for scalar field color maps
- Use `labels()` to show field values in cells

---

## Additional Resources

- Basilisk Documentation: http://basilisk.fr/

---
