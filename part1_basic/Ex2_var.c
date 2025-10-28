
/**
 * PART 1: BASIC GRID OPERATIONS IN BASILISK
 * ==========================================
 * This tutorial demonstrates:
 * - Section 2: add a variable field to the grid and visualize it
 */

#include "grid/quadtree.h"  // For adaptive mesh refinement (AMR)
#include "view.h"           // For visualization

scalar temp[];              // Declare a scalar field

int main() {
  // Define the computational domain: Domain size: [-0.5, 0.5] x [-0.5, 0.5]
  size(1.0);               // Total domain size = 1.0 unit
  origin(-0.5, -0.5);      // Set origin at center
  
  // Set initial grid resolution
  init_grid(4);            // 4x4 cells

  // Define circle parameters
  double circle_x = 0.0;      // Circle center x
  double circle_y = 0.0;      // Circle center y

  // Initialize temperature field: hot center, cool outside
  foreach() {
    double r = sqrt(sq(x - circle_x) + sq(y - circle_y));  // Distance from center
    // Hot center (30°C), cool outside (20°C)
    // Using exponential decay from center
    temp[] = 20.0 + 10.0 * exp(-60.0 * sq(r));
  }

  // Calculate field statistics
  stats s = statsf(temp);
  fprintf(stderr, "Temperature statistics:\n");
  fprintf(stderr, "  Min: %g °C\n", s.min);
  fprintf(stderr, "  Max: %g °C\n", s.max);
  fprintf(stderr, "  Mean: %g °C\n", s.sum/s.volume);
  
  // Visualize the grid
  view(width=1024, height=1024);  // Set image size
  squares("temp", min=s.min, max=s.max, linear=true);  // Color map
  cells();  // Draw grid cells
  box();  // Draw domain boundary
  labels("temp", lw=0.5);
  save("Temperature_field.png");
  dump("temp_variables");
  
  return 0;
}