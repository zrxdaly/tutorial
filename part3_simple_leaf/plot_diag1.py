#!/usr/bin/env python3
"""
Plot CO2 Flux Data from Multiple Exercises

This script reads the 'diag1' files from Ex1, Ex2, Ex3, and Ex4 folders
and plots the CO2 flux over time for comparison.

File format:
    Each diag1 file contains two columns: time and CO2_flux
    Example: "300 2.88592" means at time=300, flux=2.88592

Usage:
    python plot_diag1.py
"""

import matplotlib.pyplot as plt
import numpy as np
import os

# ============================================================================
# CONFIGURATION
# ============================================================================
# Define the exercise folders to plot
exercises = ['Ex1_NW_Hleaf', 'Ex2_W_Hleaf', 'Ex3_NW_Vleaf', 'Ex4_W_Vleaf']

# Define labels for each exercise (for the legend)
labels = [
    'Ex1: No Wind, Horizontal Leaf',
    'Ex2: With Wind, Horizontal Leaf',
    'Ex3: No Wind, Vertical Leaf',
    'Ex4: With Wind, Vertical Leaf'
]

# Define colors for each line
colors = ['blue', 'red', 'green', 'orange']

# ============================================================================
# READ DATA FROM EACH EXERCISE
# ============================================================================
print("Reading CO2 flux data from diag1 files...")

data = []  # Store data for each exercise

for i, ex in enumerate(exercises):
    # Construct the full path to the diag1 file
    filepath = os.path.join(ex, 'diag1')

    # Check if the file exists
    if os.path.exists(filepath):
        # Load data: column 0 = time, column 1 = flux
        time, flux = np.loadtxt(filepath, unpack=True)
        data.append({'time': time, 'flux': flux, 'label': labels[i], 'color': colors[i]})
        print(f"  ✓ Loaded {filepath}: {len(time)} data points")
    else:
        print(f"  ✗ Warning: {filepath} not found, skipping...")
        data.append(None)

# ============================================================================
# CREATE THE PLOT
# ============================================================================
print("\nCreating plot...")

# Create a figure with larger size for better readability
plt.figure(figsize=(10, 6))

# Plot each dataset
for i, dataset in enumerate(data):
    if dataset is not None:
        plt.plot(dataset['time'], dataset['flux'],
                 label=dataset['label'],
                 color=dataset['color'],
                 linewidth=2,
                 marker='o',
                 markersize=3,
                 markevery=10)  # Show marker every 10 points

# ============================================================================
# CUSTOMIZE THE PLOT
# ============================================================================
plt.xlabel('Time (s)', fontsize=12, fontweight='bold')
plt.ylabel('CO2 Flux (mmol/s)?', fontsize=12, fontweight='bold')
plt.title('CO2 Diffusive Flux at Leaf Surface', fontsize=14, fontweight='bold')
plt.legend(loc='best', frameon=True, shadow=True, fontsize=10)
plt.grid(True, alpha=0.3, linestyle='--')

# Add tight layout to prevent label cutoff
plt.tight_layout()

# ============================================================================
# SAVE AND DISPLAY THE PLOT
# ============================================================================
# Save the figure as a high-resolution PNG
output_file = 'co2_flux_comparison.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\n✓ Plot saved as: {output_file}")

# Display the plot
plt.show()
print("\n✓ Done!")
