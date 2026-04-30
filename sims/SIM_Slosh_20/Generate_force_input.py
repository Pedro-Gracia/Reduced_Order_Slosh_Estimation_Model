import numpy as np
import matplotlib.pyplot as plt

'''============================================================'''
'''--------------- Smooth acceleration function ---------------'''
'''============================================================'''

def a_smooth(t):
    # Transition period is between 0.5 and 3.5
    if t < 0.5:
        return 0.0
    elif t <= 3.5:
        # Intuition: Using a (1 - cos) function creates a "S-curve" 
        # that starts and ends with zero derivative (smooth ramp-up/down)
        return 0.75 * 0.5 * (1 - np.cos(np.pi * (t - 0.5) / (3.5 - 0.5)))
    else:
        return 0.0

'''============================================================'''
'''--------------- Generate Data and Table --------------------'''
'''============================================================'''

# Define the full time range for the plot
t_full = np.linspace(0.0, 5.0, 1000)
a_values = [a_smooth(ti) for ti in t_full]

# Print the table for OpenFOAM 
t_dense = np.arange(0.5, 3.5, 0.001)

print("// OpenFOAM Table Data")
print("(0.00   0.00)")
print("(0.40   0.00)")
for ti in t_dense:
    ai = a_smooth(ti)
    print(f"({ti:.3f}   {ai:.6f})")
print("(3.50   0.00)")
print("(5.00   0.00)")

'''============================================================'''
'''---------------------- Plotting ----------------------------'''
'''============================================================'''

plt.figure(figsize=(10, 6))
plt.plot(t_full, a_values, label='Smooth Acceleration $a(t)$', color='blue', lw=2)

# Formatting the plot
plt.title('Smooth Acceleration Profile for OpenFOAM')
plt.xlabel('Time (s)')
plt.ylabel('Acceleration ($m/s^2$)')
plt.grid(True, linestyle='--', alpha=0.7)
plt.axhline(0, color='black', lw=1)
plt.legend()

# Display the plot
plt.show()
