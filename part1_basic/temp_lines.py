"""
Plot temperature profile along y=0 line from Ex2_var.c

The temperature field is defined as:
temp = 20.0 + 10.0 * exp(-60.0 * r^2)
where r = sqrt((x - circle_x)^2 + (y - circle_y)^2)

For y = 0 (and circle_y = 0, circle_x = 0):
temp(x, 0) = 20.0 + 10.0 * exp(-60.0 * x^2)
"""

import numpy as np
import matplotlib.pyplot as plt

# Define the domain (from the C code)
x_min, x_max = -0.5, 0.5

# Circle parameters from the C code
circle_x = 0.0
circle_y = 0.0

# Create x values along the line y = 0
x = np.linspace(x_min, x_max, 1000)
y = 0.0  # Along y = 0 line

# Calculate temperature according to the formula
r = np.sqrt((x - circle_x)**2 + (y - circle_y)**2)
temp = 20.0 + 10.0 * np.exp(-60.0 * r**2)

# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(x, temp, 'b-', linewidth=6, label='Temperature at y=0')
# plt.axhline(y=20.0, color='r', linestyle='--', linewidth=1, alpha=0.5, label='Baseline (20°C)')
# plt.axhline(y=30.0, color='orange', linestyle='--', linewidth=1, alpha=0.5, label='Maximum (30°C)')
# plt.axvline(x=0.0, color='gray', linestyle=':', linewidth=1, alpha=0.5)

# Labels and formatting
plt.xlabel('x position', fontsize=12)
plt.ylabel('Temperature (°C)', fontsize=12)
plt.grid(True, alpha=0.3)
plt.legend(fontsize=10)
# plt.axis('off')
plt.xlim(x_min, x_max)
plt.ylim(19, 31)

# Save the figure
plt.tight_layout()
plt.savefig('temp_profile_y0.png', dpi=150, bbox_inches='tight', transparent=True)

# Display the plot
plt.show()
