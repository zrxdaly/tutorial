
/**
 * PART 1: BASIC GRID OPERATIONS IN BASILISK
 * ==========================================
 * This tutorial demonstrates:
 * - Section 4: create two hot spots
 */

#include "grid/quadtree.h"  // For adaptive mesh refinement (AMR)
#include "view.h"           // For visualization

scalar temp[];              // Declare a scalar field

int main() {
  // Define the computational domain: Domain size: [-0.5, 0.5] x [-0.5, 0.5]
  size(1.0);               // Total domain size = 1.0 unit
  origin(-0.5, -0.5);      // Set origin at center
  
  // Set initial grid resolution
  init_grid(16);            // 16x16 cells

  // Define circle parameters
  int max_level = 5;
  double circle_x1 = -0.2;      // Circle center x1
  double circle_y1 = 0.0;      // Circle center y1
  double circle_radius1 = 0.3; // Circle radius1

  // here you should define the second hot spot based on README.md
  double circle_x2 = 0.2;      // fill in Circle center x2
  double circle_y2 = 0.;      // fill in Circle center y2
  double circle_radius2 = 0.3; // Circle radius2

  // how should we modify refine function?
  refine((sq(x - circle_x1) + sq(y - circle_y1) < sq(circle_radius1) ||
          sq(x - circle_x2) + sq(y - circle_y2) < sq(circle_radius2)) && level < max_level);

  // Initialize temperature field: hot center, cool outside
  foreach() {
    double r1 = sqrt(sq(x - circle_x1) + sq(y - circle_y1));  // Distance from center1
    double r2 = sqrt(sq(x - circle_x2) + sq(y - circle_y2));  // Distance from center2
    double temp1 = 35.0 * exp(-60.0 * sq(r1));  // Hot center1 (35°C)
    double temp2 = 28.0 * exp(-30.0 * sq(r2));  // Hot center2 (28°C)

    temp[] = 18.0 + temp1 + temp2;
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








