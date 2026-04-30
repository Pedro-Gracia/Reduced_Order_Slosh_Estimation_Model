"""
============================================================
FINAL ANALYTIC SLOSH MODEL:
CURVE-FIRST HEIGHT-INTERPOLATED PARAMETER SURFACES

This script builds continuous force-model parameter surfaces of the form

    theta(t,h)

where:
    t -> time
    h -> fuel height

The important idea in this version is that the time-history curves are
constructed first at each known fuel height. Then, for every fixed time
sample, the parameter value is interpolated across fuel height.

This avoids fitting the time-model parameters directly as functions of
height, which can distort the original time dynamics.

Workflow:
    1. Load analytic time-fit parameters from each fuel-height JSON file.
    2. Reconstruct each coefficient curve theta(t,h_i) at known heights.
    3. At each fixed time t_k, interpolate theta(t_k,h_i) across height.
    4. Assemble the full surface theta(t,h).
    5. Plot and save the resulting surfaces for report/analysis use.
============================================================
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import PchipInterpolator, CubicSpline, interp1d


''''============================================================'''
''''--------------------- User Settings -------------------------'''
''''============================================================'''

# Base directory of this script
# All simulation folders and output files are built relative to this path
base_dir = os.path.dirname(__file__)

# Force-model parameters that will be converted into theta(t,h) surfaces
target_vars = ["c1", "c2", "c3", "c4", "d", "b"]

# Fuel-height simulation folders to include
# Example: idx = 30 corresponds to h = 0.30 and folder SIM_Slosh_30
height_indices = list(range(20, 90, 10))

# Time grid used to reconstruct the analytic coefficient curves
t_min = 0.0
t_max = 100.0
n_time = 600

# Number of fuel-height points used in the interpolated surface
n_height = 120

# Height interpolation method
# Recommended:
#   "pchip"  -> smooth, shape-preserving, avoids bad overshoot
# Options:
#   "pchip"  -> monotone-preserving piecewise cubic Hermite interpolation
#   "cubic"  -> cubic spline, smoother but can overshoot
#   "linear" -> safest and most conservative, but less smooth
height_interp_kind = "pchip"

# Selected times used to visualize interpolation behavior in the h direction
h_slice_times = [0.0, 1.0, 2.5, 5.0, 10.0, 25.0, 50.0, 100.0]

# Exponential clipping limit used to avoid overflow during model evaluation
MAX_EXP = 100.0


''''============================================================'''
''''--------------------- Parameter Name Maps -------------------'''
''''============================================================'''

# Defines the meaning and order of parameters for each supported analytic model
param_names_map = {
    "damped_sine_exp": ["c_inf", "A", "tau1", "omega", "phi", "B", "tau2"],
    "exp2": ["p0", "p1", "p2", "p3", "p4"]
}


''''============================================================'''
''''--------------------- Load Fit Summary Files ----------------'''
''''============================================================'''

# Storage dictionary for the raw JSON fit data
# Each parameter stores:
#   heights      -> known fuel heights
#   params       -> analytic time-model parameters at each height
#   model        -> model type used for that parameter
#   expected_len -> number of model parameters expected
raw_data = {
    var: {
        "heights": [],
        "params": [],
        "model": None,
        "expected_len": None
    }
    for var in target_vars
}

print("============================================================")
print("LOADING FIT SUMMARY JSON FILES")
print("============================================================")

# Loop through all requested simulation heights and load the JSON summaries
for idx in height_indices:
    h = idx / 100.0

    filename = f"sims/SIM_Slosh_{idx}/Analytic_Force_{idx}.json"
    filepath = os.path.join(base_dir, filename)

    # Skip missing simulation folders/files without stopping the full script
    if not os.path.exists(filepath):
        print(f"WARNING: {filename} not found. Skipping h={h:.2f}.")
        continue

    # Load analytic fit summary for the current fuel height
    with open(filepath, "r") as f:
        fit_summary = json.load(f)

    print(f"\n--- Reading h={h:.2f} ---")

    # Extract each force-model parameter from the current JSON file
    for var in target_vars:
        if var not in fit_summary:
            print(f"  WARNING: {var} not found in {filename}.")
            continue

        # Current analytic model and fitted parameter vector
        current_model = fit_summary[var]["best_model"]
        current_params = np.array(fit_summary[var]["params"], dtype=float)

        # Initialize the required model type and parameter length
        # This keeps all heights consistent for a given variable
        if raw_data[var]["model"] is None:
            raw_data[var]["model"] = current_model
            raw_data[var]["expected_len"] = len(current_params)

        # Accept data only if the model type and parameter length match
        if (
            current_model == raw_data[var]["model"]
            and len(current_params) == raw_data[var]["expected_len"]
        ):
            raw_data[var]["heights"].append(h)
            raw_data[var]["params"].append(current_params)

            print(
                f"  {var}: accepted model={current_model}, "
                f"n_params={len(current_params)}"
            )
        else:
            print(
                f"  WARNING: Skipping {var} at h={h:.2f}. "
                f"Found model={current_model}, n={len(current_params)}, "
                f"expected model={raw_data[var]['model']}, "
                f"n={raw_data[var]['expected_len']}."
            )

print("============================================================")


''''============================================================'''
''''--------------------- Build Model Structures ----------------'''
''''============================================================'''

# Main dictionary used by the rest of the script
# Each variable stores known-height data, reconstructed curves, and surfaces
models = {}

print("\n============================================================")
print("BUILDING CURVE-FIRST MODEL STRUCTURES")
print("============================================================")

for var in target_vars:

    # Convert accepted raw data into NumPy arrays
    h_raw = np.array(raw_data[var]["heights"], dtype=float)
    p_raw = np.array(raw_data[var]["params"], dtype=float)
    model_type = raw_data[var]["model"]

    # Skip variables with no valid height data
    if len(h_raw) == 0:
        print(f"Skipping {var}: no valid data.")
        continue

    # Skip unsupported analytic model types
    if model_type not in param_names_map:
        print(f"Skipping {var}: unsupported model type {model_type}.")
        continue

    # Sort all data by increasing fuel height
    sort_idx = np.argsort(h_raw)
    h_raw = h_raw[sort_idx]
    p_raw = p_raw[sort_idx, :]

    # Check that the parameter vector length matches the expected model
    expected_n = len(param_names_map[model_type])

    if p_raw.shape[1] != expected_n:
        print(
            f"Skipping {var}: parameter count mismatch. "
            f"Found {p_raw.shape[1]}, expected {expected_n}."
        )
        continue

    # Store the structured model information
    models[var] = {
        "model_type": model_type,
        "param_names": param_names_map[model_type],
        "h_raw": h_raw,
        "p_raw": p_raw,
        "curves_raw": None,
        "surface": None
    }

    print(
        f"{var}: model={model_type}, "
        f"heights={h_raw}, p_shape={p_raw.shape}"
    )

print("============================================================")


''''============================================================'''
''''--------------------- Time Model Evaluation -----------------'''
''''============================================================'''

def evaluate_time_model(t, model_type, params):
    """
    Evaluate theta(t) at one known fuel height.

    Inputs:
        t          -> time array where the coefficient is evaluated
        model_type -> analytic model form used for this coefficient
        params     -> fitted model parameters at one fuel height

    Output:
        y -> reconstructed coefficient time history theta(t,h_i)
    """

    # Make sure time is treated as a floating-point NumPy array
    t = np.asarray(t, dtype=float)

    # ------------------------------------------------------------
    # Damped sinusoid plus exponential drift model
    # ------------------------------------------------------------
    if model_type == "damped_sine_exp":
        c_inf, A, tau1, omega, phi, B, tau2 = params

        # Protect against division by zero in exponential time constants
        tau1 = tau1 if abs(tau1) >= 1.0e-8 else 1.0e-8
        tau2 = tau2 if abs(tau2) >= 1.0e-8 else 1.0e-8

        # Clip exponential arguments to avoid numerical overflow
        exp1_arg = np.clip(-t / tau1, -np.inf, MAX_EXP)
        exp2_arg = np.clip(-t / tau2, -np.inf, MAX_EXP)

        # Reconstruct the coefficient time history
        y = (
            c_inf
            + A * np.exp(exp1_arg) * np.sin(omega * t + phi)
            + B * np.exp(exp2_arg)
        )

    # ------------------------------------------------------------
    # Two-exponential model
    # ------------------------------------------------------------
    elif model_type == "exp2":
        p0, p1, p2, p3, p4 = params

        # Clip exponential arguments to avoid numerical overflow
        exp1_arg = np.clip(p1 * t, -np.inf, MAX_EXP)
        exp2_arg = np.clip(p3 * t, -np.inf, MAX_EXP)

        # Reconstruct the coefficient time history
        y = p0 * np.exp(exp1_arg) + p2 * np.exp(exp2_arg) + p4

    # ------------------------------------------------------------
    # Unsupported model type
    # ------------------------------------------------------------
    else:
        raise ValueError(f"Unsupported model type: {model_type}")

    # Replace NaN and infinite values with finite numbers for plotting/saving
    return np.nan_to_num(
        y,
        nan=0.0,
        posinf=np.finfo(float).max,
        neginf=-np.finfo(float).max
    )


''''============================================================'''
''''---------------- Interpolation Across Fuel Height -----------'''
''''============================================================'''

def interpolate_across_height(h_raw, y_raw, h_grid, kind="pchip"):
    """
    Interpolate theta(t_k,h_i) across fuel height for one fixed time t_k.

    Inputs:
        h_raw  -> known simulated fuel heights
        y_raw  -> coefficient values at those known heights for one time
        h_grid -> fuel-height grid where the coefficient is interpolated
        kind   -> interpolation method in the fuel-height direction

    Output:
        y_grid -> interpolated coefficient values across h_grid
    """

    # Convert all inputs to NumPy arrays for interpolation
    h_raw = np.asarray(h_raw, dtype=float)
    y_raw = np.asarray(y_raw, dtype=float)
    h_grid = np.asarray(h_grid, dtype=float)

    # Shape-preserving cubic Hermite interpolation
    if kind == "pchip":
        interpolator = PchipInterpolator(h_raw, y_raw, extrapolate=False)
        y_grid = interpolator(h_grid)

    # Natural cubic spline interpolation
    elif kind == "cubic":
        interpolator = CubicSpline(h_raw, y_raw, bc_type="natural", extrapolate=False)
        y_grid = interpolator(h_grid)

    # Linear interpolation
    elif kind == "linear":
        interpolator = interp1d(
            h_raw,
            y_raw,
            kind="linear",
            bounds_error=False,
            fill_value=np.nan
        )
        y_grid = interpolator(h_grid)

    # Invalid interpolation option
    else:
        raise ValueError(f"Unknown height interpolation kind: {kind}")

    return y_grid


''''============================================================'''
''''---------------- Build Curve-First Surfaces -----------------'''
''''============================================================'''

# Common time grid used to reconstruct all known-height curves
t_grid = np.linspace(t_min, t_max, n_time)

print("\n============================================================")
print("RECONSTRUCTING KNOWN-HEIGHT TIME CURVES")
print("AND INTERPOLATING ACROSS HEIGHT")
print("============================================================")

for var in target_vars:
    if var not in models:
        continue

    # Extract model information for the current parameter
    model = models[var]
    model_type = model["model_type"]
    h_raw = model["h_raw"]
    p_raw = model["p_raw"]

    # Fuel-height grid for the interpolated surface
    h_grid = np.linspace(h_raw.min(), h_raw.max(), n_height)

    # curves_raw shape:
    #   rows    -> known fuel heights
    #   columns -> time samples
    curves_raw = np.zeros((len(h_raw), len(t_grid)))

    # Reconstruct theta(t,h_i) at each known simulated height
    for i, h_i in enumerate(h_raw):
        curves_raw[i, :] = evaluate_time_model(
            t=t_grid,
            model_type=model_type,
            params=p_raw[i, :]
        )

    # Surface shape expected by plot_surface:
    #   rows    -> interpolated fuel-height grid
    #   columns -> time samples
    Z = np.zeros((len(h_grid), len(t_grid)))

    # For each fixed time, interpolate coefficient values across fuel height
    for j in range(len(t_grid)):
        y_at_known_heights = curves_raw[:, j]

        Z[:, j] = interpolate_across_height(
            h_raw=h_raw,
            y_raw=y_at_known_heights,
            h_grid=h_grid,
            kind=height_interp_kind
        )

    # Store reconstructed curves and interpolated surface
    models[var]["t_grid"] = t_grid
    models[var]["h_grid"] = h_grid
    models[var]["curves_raw"] = curves_raw
    models[var]["surface"] = Z

    print(
        f"{var}: curves_raw shape={curves_raw.shape}, "
        f"surface shape={Z.shape}, interpolation={height_interp_kind}"
    )

print("============================================================")


''''============================================================'''
''''--------------- Known-Height Surface Evaluation ------------'''
''''============================================================'''

def evaluate_surface_at_known_height(var_name, h_value):
    """
    Evaluate the curve-first surface at a requested fuel height.

    For a height that already exists in the simulated data, this should match
    the original reconstructed known-height curve very closely.
    """

    # Extract model data for the requested parameter
    model = models[var_name]

    h_raw = model["h_raw"]
    curves_raw = model["curves_raw"]

    # Initialize evaluated surface curve
    y = np.zeros(len(model["t_grid"]))

    # At each time, interpolate across height to the requested h_value
    for j in range(len(model["t_grid"])):
        y[j] = interpolate_across_height(
            h_raw=h_raw,
            y_raw=curves_raw[:, j],
            h_grid=np.array([h_value]),
            kind=height_interp_kind
        )[0]

    return y


''''============================================================'''
''''--------------- Known-Height Preservation Check ------------'''
''''============================================================'''

print("\n============================================================")
print("KNOWN-HEIGHT PRESERVATION CHECK")
print("============================================================")

# Compare the surface-evaluated curve against the original curve at each
# simulated height. For exact known heights, the error should be near zero.
for var in target_vars:
    if var not in models:
        continue

    model = models[var]

    print(f"\n--- {var} ---")

    for i, h_i in enumerate(model["h_raw"]):
        y_interp = evaluate_surface_at_known_height(var, h_i)
        y_true = model["curves_raw"][i, :]

        err = np.sqrt(np.mean((y_interp - y_true) ** 2))

        print(f"  h={h_i:.2f}: RMSE(surface vs original curve) = {err:.6e}")

print("============================================================")


''''============================================================'''
''''--------------- Print Compact Model Summary ----------------'''
''''============================================================'''

print("\n============================================================")
print("MODEL SUMMARY")
print("============================================================")

print(
    "Surface construction:\n"
    "1. Reconstruct theta(t,h_i) at each known fuel height using the time-fit JSON parameters.\n"
    "2. At each fixed time t_k, interpolate theta(t_k,h_i) across fuel height h.\n"
    "3. The final surface theta(t,h) is therefore exact at known heights and smooth between heights.\n"
)

# Print model type, available heights, and interpolation method for each parameter
for var in target_vars:
    if var not in models:
        continue

    model = models[var]

    print(f"\n{var}(t,h)")
    print(f"  time model: {model['model_type']}")
    print(f"  known heights: {model['h_raw']}")
    print(f"  height interpolation: {height_interp_kind}")

print("============================================================")


''''============================================================'''
''''--------------- Plot H-Direction Fits -----------------------'''
''''============================================================'''

def plot_h_direction_fits_for_var(var_name):
    """
    Plot the interpolation behavior in the fuel-height direction.

    For selected time slices, this function plots:
        1. known-height values theta(t_k,h_i)
        2. interpolated curve theta(t_k,h)

    This directly shows how the surface is constructed between simulated
    fuel heights.
    """

    # Extract stored model data
    model = models[var_name]

    h_raw = model["h_raw"]
    h_grid = model["h_grid"]
    t_grid_local = model["t_grid"]
    curves_raw = model["curves_raw"]

    # Set subplot layout based on number of requested time slices
    n_times = len(h_slice_times)
    ncols = 2
    nrows = int(np.ceil(n_times / ncols))

    fig, axs = plt.subplots(nrows, ncols, figsize=(14, 3.5 * nrows))
    axs = np.array(axs).reshape(-1)

    # Plot one h-direction interpolation curve for each selected time
    for k, t_request in enumerate(h_slice_times):
        ax = axs[k]

        # Find nearest available time sample on the reconstruction grid
        j = int(np.argmin(np.abs(t_grid_local - t_request)))
        t_actual = t_grid_local[j]

        # Values at known heights for this fixed time
        y_raw = curves_raw[:, j]

        # Interpolated h-direction curve for this fixed time
        y_interp = interpolate_across_height(
            h_raw=h_raw,
            y_raw=y_raw,
            h_grid=h_grid,
            kind=height_interp_kind
        )

        # Plot known-height coefficient values
        ax.plot(
            h_raw,
            y_raw,
            "o",
            markersize=7,
            label="known-height curves"
        )

        # Plot interpolated coefficient values across height
        ax.plot(
            h_grid,
            y_interp,
            "-",
            linewidth=2.2,
            label=f"Interpolation at t={t_actual:.2f} s"
        )

        # Axis labels and plot formatting
        # ax.set_title(f"{var_name}: h-fit ")
        ax.set_xlabel("Fuel height h [-]")
        ax.set_ylabel(f"{var_name}(t,h)")
        ax.grid(True, alpha=0.5)
        ax.legend(fontsize=8)

    # Remove unused subplot panels
    for j in range(n_times, len(axs)):
        fig.delaxes(axs[j])

    # Overall figure title
    fig.suptitle(
        f"Fuel height Direction Fits for {var_name} at Selected Times",
        fontsize=16
    )

    # plt.tight_layout()
    plt.subplots_adjust(top=0.92)


def plot_all_h_direction_fits():
    """
    Generate h-direction interpolation plots for every available parameter.
    """

    print("\n============================================================")
    print("GENERATING H-DIRECTION FIT PLOTS")
    print("============================================================")

    for var in target_vars:
        if var in models:
            print(f"-> Generating h-direction fits for {var}...")
            plot_h_direction_fits_for_var(var)


''''============================================================'''
''''----------- Plot Original Curves vs Surface Curves ----------'''
''''============================================================'''

def plot_known_height_curve_check_for_var(var_name):
    """
    Plot original known-height curves against surface-evaluated curves.

    This confirms that the curve-first interpolation method preserves the
    original analytic time curves at the simulated fuel heights.
    """

    # Extract stored model data
    model = models[var_name]

    h_raw = model["h_raw"]
    t_grid_local = model["t_grid"]
    curves_raw = model["curves_raw"]

    fig, ax = plt.subplots(figsize=(12, 6))

    # Assign one color per known fuel height
    colors = plt.cm.jet(np.linspace(0, 1, len(h_raw)))

    # Plot original curve and surface-evaluated curve for each known height
    for i, h_i in enumerate(h_raw):
        y_original = curves_raw[i, :]
        y_surface = evaluate_surface_at_known_height(var_name, h_i)

        ax.plot(
            t_grid_local,
            y_original,
            color=colors[i],
            linewidth=2.2,
            label=f"original h={h_i:.2f}"
        )

        ax.plot(
            t_grid_local,
            y_surface,
            "--",
            color=colors[i],
            linewidth=1.5,
            label=f"surface h={h_i:.2f}"
        )

    # Axis labels and plot formatting
    ax.set_title(f"Known-Height Preservation Check: {var_name}")
    ax.set_xlabel("Time t [s]")
    ax.set_ylabel(f"{var_name}(t,h)")
    ax.grid(True, alpha=0.5)
    ax.legend(ncol=2, fontsize=8)

    plt.tight_layout()


def plot_all_known_height_checks():
    """
    Generate known-height preservation check plots for all parameters.
    """

    print("\n============================================================")
    print("GENERATING KNOWN-HEIGHT CURVE CHECKS")
    print("============================================================")

    for var in target_vars:
        if var in models:
            print(f"-> Generating known-height check for {var}...")
            plot_known_height_curve_check_for_var(var)


''''============================================================'''
''''---------------- Plot 3D Parameter Surface ------------------'''
''''============================================================'''

def plot_3d_for_var(var_name):
    """
    Plot the full curve-first height-interpolated surface theta(t,h).

    The translucent surface shows the interpolated parameter surface, while
    the bold curves show the original known-height analytic curves used to
    construct the surface.
    """

    # Extract stored model data
    model = models[var_name]

    t_grid_local = model["t_grid"]
    h_grid = model["h_grid"]
    h_raw = model["h_raw"]
    curves_raw = model["curves_raw"]
    Z = model["surface"]

    # Build meshgrid for 3D plotting
    T, H = np.meshgrid(t_grid_local, h_grid)

    # Create 3D figure and axis
    fig = plt.figure(figsize=(11, 8))
    ax = fig.add_subplot(111, projection="3d")

    # Plot interpolated theta(t,h) surface
    ax.plot_surface(
        T,
        H,
        Z,
        cmap=cm.plasma,
        edgecolor="none",
        alpha=0.50
    )

    # Assign one color per known fuel height
    colors = plt.cm.jet(np.linspace(0, 1, len(h_raw)))

    # Determine z-axis limits from the surface values
    z_min = np.nanmin(Z)
    z_max = np.nanmax(Z)

    # Overlay the original known-height curves on top of the surface
    for i, h_i in enumerate(h_raw):
        h_line = np.full_like(t_grid_local, h_i)

        ax.plot(
            t_grid_local,
            h_line,
            curves_raw[i, :],
            color=colors[i],
            linewidth=4.0,
            zorder=10,
            label=f"h={h_i:.2f}"
        )

        ax.scatter(
            t_grid_local[::50],
            h_line[::50],
            curves_raw[i, ::50],
            color=colors[i],
            s=18,
            zorder=11
        )

    # Add a small z-axis margin so the curves are visible
    z_margin = 0.10 * (z_max - z_min)

    if not np.isfinite(z_margin) or z_margin == 0.0:
        z_margin = 1.0

    ax.set_zlim(z_min - z_margin, z_max + z_margin)

    # Axis labels and title
    ax.set_xlabel("Time t [s]")
    ax.set_ylabel("Fuel height h [-]")
    ax.set_zlabel(f"{var_name}(t,h)")
    ax.set_title(f"Curve-First Height-Interpolated Surface: {var_name}(t,h)")

    # Set 3D plot proportions and viewing angle
    ax.set_box_aspect([1.6, 1.0, 0.8])
    ax.view_init(elev=30, azim=-135)
    ax.legend()

    # Add compact equation annotation to the plot
    equation_text = (
        rf"${var_name}(t,h):\ "
        r"\theta(t,h_i)\ \mathrm{from\ time\ fits},\quad "
        r"\theta(t,h)=\mathcal{I}_h\{\theta(t,h_i)\}$"
    )

    ax.text2D(
        0.03,
        0.96,
        equation_text,
        transform=ax.transAxes,
        fontsize=11,
        bbox=dict(
            boxstyle="round,pad=0.5",
            facecolor="white",
            edgecolor="gray",
            alpha=0.85
        )
    )

    plt.tight_layout()


def plot_all_3d_surfaces():
    """
    Generate 3D theta(t,h) surfaces for all available parameters.
    """

    print("\n============================================================")
    print("GENERATING 3D SURFACES")
    print("============================================================")

    for var in target_vars:
        if var in models:
            print(f"-> Generating 3D surface for {var}...")
            plot_3d_for_var(var)


''''============================================================'''
''''---------------- Save Surface Data --------------------------'''
''''============================================================'''

# Output file containing all reconstructed curves and interpolated surfaces
outfile_npz = os.path.join(base_dir, "curve_first_height_interpolated_surfaces.npz")

# Dictionary that will be saved into the NPZ output file
save_data = {}

# Save all available surface data for each parameter
for var in target_vars:
    if var not in models:
        continue

    save_data[f"{var}_t_grid"] = models[var]["t_grid"]
    save_data[f"{var}_h_grid"] = models[var]["h_grid"]
    save_data[f"{var}_h_raw"] = models[var]["h_raw"]
    save_data[f"{var}_curves_raw"] = models[var]["curves_raw"]
    save_data[f"{var}_surface"] = models[var]["surface"]
    save_data[f"{var}_params_raw"] = models[var]["p_raw"]

# Write the curve-first height-interpolated surface data to disk
np.savez(outfile_npz, **save_data)

print(f"\nSaved curve-first interpolated surfaces to:\n{outfile_npz}")



''''============================================================'''
''''---------------- Print LaTeX Model Equations ----------------'''
''''============================================================'''

def print_report_equations():
    """
    Print LaTeX-ready analytic model equations for the report.

    For each force-model coefficient, this prints:
        1. the analytic time model at each known fuel height
        2. the final interpolated model theta(t,h)
        3. the interpolation definition used across height
    """

    print("\n============================================================")
    print("LATEX-READY MODEL EQUATIONS FOR REPORT")
    print("============================================================\n")

    for var in target_vars:
        if var not in models:
            continue

        model = models[var]
        h_raw = model["h_raw"]
        params = model["p_raw"]
        model_type = model["model_type"]

        print(f"\n% ==========================================================")
        print(f"% {var.upper()} MODEL")
        print(f"% ==========================================================\n")

        # ----------------------------------------------------
        # Time models at each height
        # ----------------------------------------------------
        print(f"% --- Time-domain models at each height ---\n")

        for i, h_i in enumerate(h_raw):
            p = params[i]

            # Print damped-sine plus exponential-drift model
            if model_type == "damped_sine_exp":
                c_inf, A, tau1, omega, phi, B, tau2 = p

                print(r"\begin{equation}")
                print(
                    rf"{var}(t, h={h_i:.2f}) = "
                    rf"{c_inf:.6f} "
                    rf"+ ({A:.6f}) e^{{-t/{tau1:.6f}}}"
                    rf"\sin({omega:.6f}t + {phi:.6f}) "
                    rf"+ ({B:.6f}) e^{{-t/{tau2:.6f}}}"
                )
                print(r"\end{equation}\n")

            # Print two-exponential model
            elif model_type == "exp2":
                p0, p1, p2, p3, p4 = p

                print(r"\begin{equation}")
                print(
                    rf"{var}(t, h={h_i:.2f}) = "
                    rf"{p0:.6f} e^{{{p1:.6f} t}} "
                    rf"+ {p2:.6f} e^{{{p3:.6f} t}} "
                    rf"+ {p4:.6f}"
                )
                print(r"\end{equation}\n")

        # ----------------------------------------------------
        # Final interpolated model
        # ----------------------------------------------------
        print(f"% --- Final interpolated model ---\n")

        print(r"\begin{equation}")
        print(
            rf"{var}(t,h) = \mathcal{{I}}_h\left\{{ {var}(t,h_i) \right\}}"
        )
        print(r"\end{equation}")

        print("\nWhere:")
        print(r"\begin{itemize}")
        print(r"\item $h_i \in \{ " + ", ".join([f"{h:.2f}" for h in h_raw]) + r" \}$")
        print(r"\item $\mathcal{I}_h$ denotes interpolation in fuel height")
        print(r"\item Interpolation method: " + height_interp_kind)
        print(r"\end{itemize}\n")

    print("============================================================")


''''============================================================'''
''''--------------------- Main Execution ------------------------'''
''''============================================================'''

print("\n============================================================")
print("RUNNING CURVE-FIRST HEIGHT-INTERPOLATED MODEL")
print("============================================================")

# Print report-ready equations
print_report_equations()

# Generate plots showing interpolation across fuel height at fixed times
plot_all_h_direction_fits()

# Generate plots checking that the surface preserves known-height curves
plot_all_known_height_checks()

# Generate 3D parameter surfaces
plot_all_3d_surfaces()

# Display all generated figures
plt.show()
