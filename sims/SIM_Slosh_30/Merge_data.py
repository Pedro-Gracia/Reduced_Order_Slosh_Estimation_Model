import numpy as np
import os
from glob import glob
import matplotlib.pyplot as plt

''''============================================================'''
''''--------------- Set Base Case File Paths --------------------'''
''''============================================================'''

# Get the absolute path of the current Python script
# This makes the script portable because all other paths are built relative to it
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Path to the OpenFOAM force output folders
FORCES_PATH = os.path.join(BASE_DIR, "postProcessing", "forces")

# Path to the liquid center-of-mass output file
COM_PATH    = os.path.join(BASE_DIR, "postProcessing", "liquidCOM", "0", "liquidCOM.dat")


''''============================================================'''
''''--------------- Force File Reader Function ------------------'''
''''============================================================'''

def read_forces(file_path):
    """
    Read one OpenFOAM forces.dat file and extract the x-direction force.

    The OpenFOAM forces file contains pressure and viscous force components.
    For this project, the total wall force in the x-direction is computed as:

        Fx = Fx_pressure + Fx_viscous

    Inputs:
        file_path -> path to one forces.dat file

    Outputs:
        t  -> time samples from the force file
        Fx -> total x-direction force at each time sample
    """

    # Initialize storage lists for time and force
    t, Fx = [], []

    # Open and read the OpenFOAM force file line-by-line
    with open(file_path, 'r') as f:
        for line in f:

            # Skip header lines and blank lines
            if line.startswith('#') or len(line.strip()) == 0:
                continue

            # Remove parentheses from OpenFOAM vector formatting
            clean = line.replace('(', '').replace(')', '')

            # Split the cleaned line into individual numeric entries
            vals = clean.split()

            # Skip incomplete lines for safety
            if len(vals) < 7:
                continue

            # Extract simulation time
            time = float(vals[0])

            # Extract total x-force as pressure contribution plus viscous contribution
            fx = float(vals[1]) + float(vals[4])  # pressure + viscous

            # Store current time and force value
            t.append(time)
            Fx.append(fx)

    return np.array(t), np.array(Fx)


''''============================================================'''
''''--------------- Load and Combine Force Data -----------------'''
''''============================================================'''

# Find all available OpenFOAM forces.dat files
# This supports cases where OpenFOAM writes force data into multiple time folders
force_files = sorted(glob(os.path.join(FORCES_PATH, "*", "forces.dat")))

# Stop the script if no force files are found
if len(force_files) == 0:
    print("DEBUG: searched path:", FORCES_PATH)
    raise RuntimeError("No forces.dat found")

# Print all detected force files for verification
print("\nFound force files:")
for f in force_files:
    print("  ", f)

# Initialize storage for all force-file segments
t_all = []
Fx_all = []

# Read each force file and append its valid data
for f in force_files:
    t, Fx = read_forces(f)

    if len(t) > 0:
        t_all.append(t)
        Fx_all.append(Fx)

# Combine all force segments into one continuous force history
t_force = np.concatenate(t_all)
Fx_force = np.concatenate(Fx_all)

# Sort force data by time
# This is important when multiple postProcessing folders are present
idx = np.argsort(t_force)
t_force = t_force[idx]
Fx_force = Fx_force[idx]

# Remove duplicate force samples
# This prevents repeated times from causing alignment problems later
mask = np.diff(t_force, prepend=-1) > 1e-10
t_force = t_force[mask]
Fx_force = Fx_force[mask]

# Print total number of valid force samples
print(f"\nTotal force samples: {len(t_force)}")


''''============================================================'''
''''--------------- Load Liquid COM Data ------------------------'''
''''============================================================'''

# Stop the script if the liquid COM file does not exist
if not os.path.exists(COM_PATH):
    raise RuntimeError("COM file not found")

# Print the COM file being read
print("\nReading COM:", COM_PATH)

# Load liquid center-of-mass output from OpenFOAM
com_data = np.loadtxt(COM_PATH, comments="#")

# Extract COM time and x-position histories
t_com = com_data[:, 0]
x_com = com_data[:, 1]

# Sort COM data by time for safety
idx = np.argsort(t_com)
t_com = t_com[idx]
x_com = x_com[idx]

# Remove duplicate COM samples
mask = np.diff(t_com, prepend=-1) > 1e-10
t_com = t_com[mask]
x_com = x_com[mask]

# Print total number of valid COM samples
print(f"Total COM samples: {len(t_com)}")


''''============================================================'''
''''--------------- Align Force and COM Data --------------------'''
''''============================================================'''

# Find time samples that exist in both the force data and COM data
# The times are rounded to avoid tiny floating-point mismatch issues
t_common = np.intersect1d(
    np.round(t_force, 6),
    np.round(t_com, 6)
)

# Initialize aligned force and COM arrays
Fx_aligned = []
x_aligned = []

# Loop through common time samples and extract matching force and COM values
for t_val in t_common:
    idx_f = np.where(np.isclose(t_force, t_val))[0][0]
    idx_c = np.where(np.isclose(t_com, t_val))[0][0]

    Fx_aligned.append(Fx_force[idx_f])
    x_aligned.append(x_com[idx_c])

# Convert aligned lists into NumPy arrays
t_common = np.array(t_common)
Fx_aligned = np.array(Fx_aligned)
x_aligned = np.array(x_aligned)

# Print number of samples available after alignment
print(f"Aligned samples: {len(t_common)}")


''''============================================================'''
''''--------------- Compute Relative COM Displacement -----------'''
''''============================================================'''

# Use the pre-forcing interval to estimate the equilibrium COM location
mask_eq = t_common < 0.5   # before forcing starts

# Initial/equilibrium COM position
x0 = np.mean(x_aligned[mask_eq])

# Relative COM displacement used by the estimator
qx_aligned = x_aligned - x0


''''============================================================'''
''''--------------- Plot Raw Force and COM Data -----------------'''
''''============================================================'''

fig1, ax1 = plt.subplots(2, 1, figsize=(12, 10))

# Plot raw wall force from OpenFOAM
ax1[0].plot(t_force, Fx_force, '.', label="Fx raw", color='tab:blue')
ax1[0].set_title("Raw Data: Force ($F_x$)")
ax1[0].set_ylabel("Force")
ax1[0].legend(loc="upper right")
ax1[0].grid(True)

# Plot raw liquid center-of-mass x-position from OpenFOAM
ax1[1].plot(t_com, x_com, '.', label="xCOM raw", color='tab:orange')
ax1[1].set_title("Raw Data: Center of Mass ($x_{COM}$)")
ax1[1].set_xlabel("Time [s]")
ax1[1].set_ylabel("Position")
ax1[1].legend(loc="upper right")
ax1[1].grid(True)

# Adjust spacing between subplot labels and titles
plt.tight_layout() # Prevents overlap of labels and titles


''''============================================================'''
''''--------------- Plot Aligned Estimator Signals --------------'''
''''============================================================'''

fig2, ax2 = plt.subplots(2, 1, figsize=(12, 10))

# Plot aligned force signal at the common force/COM time samples
ax2[0].plot(t_common, Fx_aligned, label="Fx (aligned)", color='green')
ax2[0].set_title("Aligned Signals: Force")
ax2[0].set_ylabel("Force [N]")
ax2[0].legend(loc="upper right")
ax2[0].grid(True)

# Plot aligned relative COM displacement used as q measurement
ax2[1].plot(t_common, x_aligned - x0, label="q (aligned)", color='red')
ax2[1].set_title("Aligned Signals: Relative Displacement ($q$)")
ax2[1].set_xlabel("Time [s]")
ax2[1].set_ylabel("Value")
ax2[1].legend(loc="upper right")
ax2[1].grid(True)

# Adjust spacing between subplot labels and titles
plt.tight_layout()


''''============================================================'''
''''--------------- Display All Figures -------------------------'''
''''============================================================'''

# Show all generated figures
# When running as a script, this will keep the figure windows open
# until they are manually closed
plt.show()


''''============================================================'''
''''--------------- Save Merged CFD Dataset ---------------------'''
''''============================================================'''

# Save the aligned CFD-derived measurements used by the estimator
# Output columns:
#   column 0 -> common time samples
#   column 1 -> aligned wall force Fx
#   column 2 -> aligned relative COM displacement qx
np.savetxt(
    os.path.join(BASE_DIR, "CFD_data.dat"),
    np.column_stack((t_common, Fx_aligned, qx_aligned)),
    header="time Fx qx",
    comments=''
)

# Print confirmation that the estimator input file was written
print("\nSaved: CFD_data.dat")
