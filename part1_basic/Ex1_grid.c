
/**
 * PART 1: BASIC GRID OPERATIONS IN BASILISK
 * ==========================================
 * This tutorial demonstrates:
 * - Section 1: Simple uniform 2D grid
 */

#include "grid/quadtree.h"  // For adaptive mesh refinement (AMR)
#include "view.h"           // For visualization

/**
 * SECTION 1: UNIFORM GRID GENERATION
 * ===================================
 * We'll create a simple 2D domain and visualize the grid structure
 */

int main() {
  // Define the computational domain: Domain size: [-0.5, 0.5] x [-0.5, 0.5]
  size(1.0);               // Total domain size = 1.0 unit
  origin(-0.5, -0.5);      // Set origin at center
  
  // Set initial grid resolution
  init_grid(4);            // 4x4 cells
  
  // Visualize the grid
  view(width=800, height=800);  // Set image size
  cells();                 // Draw grid cells
  box();                   // Draw domain boundary
  save("Uniform_grid.png");
  
  return 0;
}