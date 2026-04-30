import numpy as np
import os
from glob import glob
import matplotlib.pyplot as plt

''''============================================================'''
''''--------------- SET BASE PATH -------------------------------'''
''''============================================================'''

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

FORCES_PATH = os.path.join(BASE_DIR, "postProcessing", "forces")
COM_PATH    = os.path.join(BASE_DIR, "postProcessing", "liquidCOM", "0", "liquidCOM.dat")

''''============================================================'''
''''--------------- READ FORCES --------------------------------'''
''''============================================================'''

def read_forces(file_path):
    t, Fx = [], []

    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#') or len(line.strip()) == 0:
                continue

            clean = line.replace('(', '').replace(')', '')
            vals = clean.split()

            if len(vals) < 7:
                continue

            time = float(vals[0])
            fx = float(vals[1]) + float(vals[4])  # pressure + viscous

            t.append(time)
            Fx.append(fx)

    return np.array(t), np.array(Fx)


''''============================================================'''
''''--------------- LOAD FORCES --------------------------------'''
''''============================================================'''

force_files = sorted(glob(os.path.join(FORCES_PATH, "*", "forces.dat")))

if len(force_files) == 0:
    print("DEBUG: searched path:", FORCES_PATH)
    raise RuntimeError("No forces.dat found")

print("\nFound force files:")
for f in force_files:
    print("  ", f)

t_all = []
Fx_all = []

for f in force_files:
    t, Fx = read_forces(f)

    if len(t) > 0:
        t_all.append(t)
        Fx_all.append(Fx)

t_force = np.concatenate(t_all)
Fx_force = np.concatenate(Fx_all)

# --- SORT + REMOVE DUPLICATES (critical)
idx = np.argsort(t_force)
t_force = t_force[idx]
Fx_force = Fx_force[idx]

mask = np.diff(t_force, prepend=-1) > 1e-10
t_force = t_force[mask]
Fx_force = Fx_force[mask]

print(f"\nTotal force samples: {len(t_force)}")


''''============================================================'''
''''--------------- LOAD COM -----------------------------------'''
''''============================================================'''

if not os.path.exists(COM_PATH):
    raise RuntimeError("COM file not found")

print("\nReading COM:", COM_PATH)

com_data = np.loadtxt(COM_PATH, comments="#")

t_com = com_data[:, 0]
x_com = com_data[:, 1]

# --- sort + unique (safety)
idx = np.argsort(t_com)
t_com = t_com[idx]
x_com = x_com[idx]

mask = np.diff(t_com, prepend=-1) > 1e-10
t_com = t_com[mask]
x_com = x_com[mask]

print(f"Total COM samples: {len(t_com)}")


''''============================================================'''
''''--------------- ALIGN DATA (NO INTERPOLATION) ---------------'''
''''============================================================'''

# match times up to tolerance
t_common = np.intersect1d(
    np.round(t_force, 6),
    np.round(t_com, 6)
)

Fx_aligned = []
x_aligned = []

for t_val in t_common:
    idx_f = np.where(np.isclose(t_force, t_val))[0][0]
    idx_c = np.where(np.isclose(t_com, t_val))[0][0]

    Fx_aligned.append(Fx_force[idx_f])
    x_aligned.append(x_com[idx_c])

t_common = np.array(t_common)
Fx_aligned = np.array(Fx_aligned)
x_aligned = np.array(x_aligned)

print(f"Aligned samples: {len(t_common)}")


''''============================================================'''
''''--------------- DISPLACEMENT -------------------------------'''
''''============================================================'''

mask_eq = t_common < 0.5   # before forcing starts
x0 = np.mean(x_aligned[mask_eq])
qx_aligned = x_aligned - x0

''''============================================================'''
''''--------------- PLOTS ---------------------------------------'''
''''============================================================'''

fig1, ax1 = plt.subplots(2, 1, figsize=(12, 10))

# Subplot 1 (Top): Raw Force
ax1[0].plot(t_force, Fx_force, '.', label="Fx raw", color='tab:blue')
ax1[0].set_title("Raw Data: Force ($F_x$)")
ax1[0].set_ylabel("Force")
ax1[0].legend(loc="upper right")
ax1[0].grid(True)

# Subplot 2 (Bottom): Raw COM
ax1[1].plot(t_com, x_com, '.', label="xCOM raw", color='tab:orange')
ax1[1].set_title("Raw Data: Center of Mass ($x_{COM}$)")
ax1[1].set_xlabel("Time [s]")
ax1[1].set_ylabel("Position")
ax1[1].legend(loc="upper right")
ax1[1].grid(True)

plt.tight_layout() # Prevents overlap of labels and titles


# ==========================================
# FIGURE 2: ALIGNED SIGNALS (Two Stacked Plots)
# ==========================================
fig2, ax2 = plt.subplots(2, 1, figsize=(12, 10))

# Subplot 1 (Top): Aligned Force
ax2[0].plot(t_common, Fx_aligned, label="Fx (aligned)", color='green')
ax2[0].set_title("Aligned Signals: Force")
ax2[0].set_ylabel("Force [N]")
ax2[0].legend(loc="upper right")
ax2[0].grid(True)

# Subplot 2 (Bottom): Aligned Displacement
ax2[1].plot(t_common, x_aligned - x0, label="q (aligned)", color='red')
ax2[1].set_title("Aligned Signals: Relative Displacement ($q$)")
ax2[1].set_xlabel("Time [s]")
ax2[1].set_ylabel("Value")
ax2[1].legend(loc="upper right")
ax2[1].grid(True)

plt.tight_layout()

# ==========================================
# FINAL RENDER
# ==========================================
# This command displays all open figures. 
# If running as a script, it will pause until windows are closed.
plt.show()

''''============================================================'''
''''--------------- SAVE (FOR ESTIMATION) -----------------------'''
''''============================================================'''

np.savetxt(
    os.path.join(BASE_DIR, "CFD_data.dat"),
    np.column_stack((t_common, Fx_aligned, qx_aligned)),
    header="time Fx qx",
    comments=''
)

print("\nSaved: CFD_data.dat")
