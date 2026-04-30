''''============================================================'''
''''--------------------- Imports -------------------------------'''
''''============================================================'''

import numpy as np
import os
import json
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit


''''============================================================'''
''''---------------- Load Estimator Result Data -----------------'''
''''============================================================'''

# Define the fuel level associated with this fitting and reconstruction run
fuel_level = 60

# Set up the directory path and build the estimator-results filename
base_dir = os.path.dirname(__file__)
filename = os.path.join(base_dir, f"Estimator_results_{fuel_level}.npz")

# Load the saved MHE estimator results
data = np.load(filename)

# Extract time history using the fuel-level dependent key
t  = data[f"t_{fuel_level}"]

# Extract estimated modal parameter histories
w1 = data[f"w1_{fuel_level}"]
z1 = data[f"z1_{fuel_level}"]
g1 = data[f"g1_{fuel_level}"]
w2 = data[f"w2_{fuel_level}"]
z2 = data[f"z2_{fuel_level}"]
g2 = data[f"g2_{fuel_level}"]

# Extract estimated force-model coefficient histories
c1 = data[f"c1_{fuel_level}"]
c2 = data[f"c2_{fuel_level}"]
c3 = data[f"c3_{fuel_level}"]
c4 = data[f"c4_{fuel_level}"]
d  = data[f"d_{fuel_level}"]
b  = data[f"b_{fuel_level}"]

# Extract reconstructed COM and force histories from the estimator output
q_meas  = data[f"q_{fuel_level}"]
Fx_meas = data[f"Fx_{fuel_level}"]

# Extract estimated modal states from the MHE output
q1_est  = data[f"q1_{fuel_level}"]
q1d_est = data[f"q1d_{fuel_level}"]
q2_est  = data[f"q2_{fuel_level}"]
q2d_est = data[f"q2d_{fuel_level}"]


''''============================================================'''
''''--------------------- Helper Functions ----------------------'''
''''============================================================'''

# Compute the root-mean-square error between two signals
def rmse(y_true, y_pred):
    return np.sqrt(np.mean((y_true - y_pred)**2))


# Compute the average value over the final fraction of a signal
# This is used to estimate the asymptotic value for the fit initialization
def tail_mean(y, frac=0.25):
    i0 = int((1.0 - frac) * len(y))
    return np.mean(y[i0:])


# Save the selected model type, RMSE, and fitted parameters to a JSON file
def save_fit_summary(filename_json, fit_summary):
    serializable = {}
    for key, val in fit_summary.items():
        serializable[key] = {
            "best_model": val["best_model"],
            "rmse": float(val["rmse"]),
            "params": [float(x) for x in val["params"]]
        }

    with open(filename_json, "w") as f:
        json.dump(serializable, f, indent=2)


''''============================================================'''
''''------------------ Input Acceleration Model -----------------'''
''''============================================================'''

# IMPORTANT:
# Keep this acceleration profile identical to the MHE estimator file.
# The analytic force reconstruction must use the same known input a_x(t)
# that was used when estimating the modal states and force coefficients.

def a_smooth(ti):
    if ti < 0.5:
        return 0.0
    elif ti <= 3.5:
        return 0.75 * 0.5 * (1 - np.cos(np.pi * (ti - 0.5) / (2.0 - 0.5)))
    else:
        return 0.0

# Evaluate the known lateral acceleration input at each saved time sample
a_x = np.array([a_smooth(ti) for ti in t])


''''============================================================'''
''''---------------- Analytic Fit Model Library -----------------'''
''''============================================================'''

# 1) Constant model
# This represents a parameter that is approximately time-invariant
def model_const(t, c0):
    return c0 + 0.0*t


# 2) Single exponential decay to an asymptotic value
# This represents a parameter with one dominant transient timescale
def model_exp1(t, c_inf, A1, tau1):
    return c_inf + A1*np.exp(-t/tau1)


# 3) Double exponential decay to an asymptotic value
# This represents a parameter with both fast and slow transient behavior
def model_exp2(t, c_inf, A1, tau1, A2, tau2):
    return c_inf + A1*np.exp(-t/tau1) + A2*np.exp(-t/tau2)


# 4) Damped sinusoid about an asymptotic value
# This represents an oscillatory coefficient response with decay
def model_damped_sine(t, c_inf, A, tau, omega, phi):
    return c_inf + A*np.exp(-t/tau)*np.sin(omega*t + phi)


# 5) Damped sinusoid plus exponential drift
# This captures both oscillatory behavior and non-oscillatory transient drift
def model_damped_sine_exp(t, c_inf, A, tau, omega, phi, B, tau2):
    return (
        c_inf
        + A*np.exp(-t/tau)*np.sin(omega*t + phi)
        + B*np.exp(-t/tau2)
    )


''''============================================================'''
''''--------------------- Fit Wrapper Logic ---------------------'''
''''============================================================'''

def try_fit(model_name, t, y):
    """
    Returns:
        params, y_fit, fit_rmse
    Raises:
        RuntimeError if fit fails
    """

    # Estimate a reasonable asymptotic value from the tail of the signal
    y_tail = tail_mean(y, frac=0.25)

    # Initial value of the signal
    y0 = y[0]

    # Signal amplitude scale used to initialize oscillatory models
    amp = max(np.max(y) - np.min(y), 1e-6)

    # ------------------------------------------------------------
    # Constant model fit
    # ------------------------------------------------------------
    if model_name == "const":
        p0 = [y_tail]
        params, _ = curve_fit(model_const, t, y, p0=p0, maxfev=20000)
        y_fit = model_const(t, *params)

    # ------------------------------------------------------------
    # Single exponential model fit
    # ------------------------------------------------------------
    elif model_name == "exp1":
        p0 = [y_tail, y0 - y_tail, 10.0]
        lb = [-np.inf, -np.inf, 1e-6]
        ub = [ np.inf,  np.inf, np.inf]
        params, _ = curve_fit(
            model_exp1, t, y, p0=p0, bounds=(lb, ub), maxfev=50000
        )
        y_fit = model_exp1(t, *params)

    # ------------------------------------------------------------
    # Double exponential model fit
    # ------------------------------------------------------------
    elif model_name == "exp2":
        p0 = [y_tail, 0.7*(y0 - y_tail), 3.0, 0.3*(y0 - y_tail), 30.0]
        lb = [-np.inf, -np.inf, 1e-6, -np.inf, 1e-6]
        ub = [ np.inf,  np.inf, np.inf,  np.inf, np.inf]
        params, _ = curve_fit(
            model_exp2, t, y, p0=p0, bounds=(lb, ub), maxfev=100000
        )
        y_fit = model_exp2(t, *params)

    # ------------------------------------------------------------
    # Damped sinusoid model fit
    # ------------------------------------------------------------
    elif model_name == "damped_sine":
        p0 = [y_tail, 0.5*amp, 8.0, 2.5, 0.0]
        lb = [-np.inf, -np.inf, 1e-6, 0.0, -2*np.pi]
        ub = [ np.inf,  np.inf, np.inf, np.inf,  2*np.pi]
        params, _ = curve_fit(
            model_damped_sine, t, y, p0=p0, bounds=(lb, ub), maxfev=100000
        )
        y_fit = model_damped_sine(t, *params)

    # ------------------------------------------------------------
    # Damped sinusoid plus exponential drift model fit
    # ------------------------------------------------------------
    elif model_name == "damped_sine_exp":
        p0 = [y_tail, 0.5*amp, 8.0, 2.5, 0.0, y0 - y_tail, 25.0]
        lb = [-np.inf, -np.inf, 1e-6, 0.0, -2*np.pi, -np.inf, 1e-6]
        ub = [ np.inf,  np.inf, np.inf, np.inf,  2*np.pi,  np.inf, np.inf]
        params, _ = curve_fit(
            model_damped_sine_exp, t, y, p0=p0, bounds=(lb, ub), maxfev=150000
        )
        y_fit = model_damped_sine_exp(t, *params)

    # ------------------------------------------------------------
    # Invalid model name
    # ------------------------------------------------------------
    else:
        raise ValueError(f"Unknown model_name = {model_name}")

    # Compute the fit error for the selected model
    fit_rmse = rmse(y, y_fit)
    return params, y_fit, fit_rmse


# Try all candidate models for one parameter and keep the model with lowest RMSE
def select_best_model(name, t, y, candidate_models):
    best = None

    for model_name in candidate_models:
        try:
            params, y_fit, fit_rmse = try_fit(model_name, t, y)

            if (best is None) or (fit_rmse < best["rmse"]):
                best = {
                    "best_model": model_name,
                    "params": params,
                    "y_fit": y_fit,
                    "rmse": fit_rmse
                }

        except Exception as e:
            print(f"[{name}] model {model_name} failed: {e}")

    if best is None:
        raise RuntimeError(f"All fits failed for parameter {name}")

    print(f"[{name}] best model = {best['best_model']}, RMSE = {best['rmse']:.6e}")
    return best


''''============================================================'''
''''---------------- Candidate Models Per Parameter -------------'''
''''============================================================'''

# Force-side parameters only.
#
# The coefficients c1-c4, d, and b are fitted using compact analytic functions.
# These analytic functions replace the raw time histories with smooth parameter
# trajectories that can later be used for force reconstruction.
#
# In this version, each force-side parameter is assigned the damped-sine plus
# exponential-drift model.

candidate_map = {
    "c1": ["damped_sine_exp"],
    "c2": ["damped_sine_exp"],
    "c3": ["damped_sine_exp"],
    "c4": ["damped_sine_exp"],
    "d" : ["damped_sine_exp"],
    "b" : ["damped_sine_exp"],
}


''''============================================================'''
''''---------------- Fit Force-Model Parameters -----------------'''
''''============================================================'''

# Dictionary used to store the selected model, fitted parameters, fitted
# trajectory, and RMSE for each force-model coefficient
fit_summary = {}

# Fit each force-model coefficient independently
fit_summary["c1"] = select_best_model("c1", t, c1, candidate_map["c1"])
fit_summary["c2"] = select_best_model("c2", t, c2, candidate_map["c2"])
fit_summary["c3"] = select_best_model("c3", t, c3, candidate_map["c3"])
fit_summary["c4"] = select_best_model("c4", t, c4, candidate_map["c4"])
fit_summary["d"]  = select_best_model("d",  t, d,  candidate_map["d"])
fit_summary["b"]  = select_best_model("b",  t, b,  candidate_map["b"])

# Extract fitted coefficient trajectories from the fit summary
c1_fit = fit_summary["c1"]["y_fit"]
c2_fit = fit_summary["c2"]["y_fit"]
c3_fit = fit_summary["c3"]["y_fit"]
c4_fit = fit_summary["c4"]["y_fit"]
d_fit  = fit_summary["d"]["y_fit"]
b_fit  = fit_summary["b"]["y_fit"]


''''============================================================'''
''''------------ Force Reconstruction Using Raw Parameters ------'''
''''============================================================'''

# Reconstruct the wall force using the original estimated coefficients
# directly from the MHE output
Fx_raw = (
    c1*q1_est
    + c2*q1d_est
    + c3*q2_est
    + c4*q2d_est
    + d*a_x
    + b
)

# Compute force reconstruction error using raw estimated parameters
rmse_Fx_raw = rmse(Fx_meas, Fx_raw)


''''============================================================'''
''''---------- Force Reconstruction Using Analytic Fits ---------'''
''''============================================================'''

# Reconstruct the wall force using the analytic fitted coefficient histories
Fx_fitfunc = (
    c1_fit*q1_est
    + c2_fit*q1d_est
    + c3_fit*q2_est
    + c4_fit*q2d_est
    + d_fit*a_x
    + b_fit
)

# Compute force reconstruction error using analytic parameter fits
rmse_Fx_fitfunc = rmse(Fx_meas, Fx_fitfunc)


''''============================================================'''
''''---------------- Save Analytic Force Results ----------------'''
''''============================================================'''

# Output files for the analytic force reconstruction
outfile_npz = os.path.join(base_dir, f"Analytic_Force_{fuel_level}.npz")
outfile_json = os.path.join(base_dir, f"Analytic_Force_{fuel_level}.json")

# Save raw and fitted coefficient trajectories, reconstructed forces,
# estimated modal states, and acceleration input
np.savez(
    outfile_npz,
    t=t,
    Fx_meas=Fx_meas,
    Fx_raw=Fx_raw,
    Fx_fitfunc=Fx_fitfunc,
    q1_est=q1_est,
    q1d_est=q1d_est,
    q2_est=q2_est,
    q2d_est=q2d_est,
    c1_raw=c1, c1_fit=c1_fit,
    c2_raw=c2, c2_fit=c2_fit,
    c3_raw=c3, c3_fit=c3_fit,
    c4_raw=c4, c4_fit=c4_fit,
    d_raw=d,   d_fit=d_fit,
    b_raw=b,   b_fit=b_fit,
    a_x=a_x
)

# Save a compact JSON summary of the selected models and fit coefficients
save_fit_summary(outfile_json, fit_summary)

# Print output locations and reconstruction errors
print("")
print(f"Saved fit trajectories to: {outfile_npz}")
print(f"Saved fit summary to:      {outfile_json}")
print("")
print(f"RMSE Fx using raw parameters      = {rmse_Fx_raw:.6e}")
print(f"RMSE Fx using analytic fits       = {rmse_Fx_fitfunc:.6e}")


''''============================================================'''
''''---------------- Plot Analytic Parameter Fits ---------------'''
''''============================================================'''

fig, axs = plt.subplots(3, 2, figsize=(14, 12), sharex=True)

# ------------------------------------------------------------
# Row 1: c1 and c2 parameter fits
# ------------------------------------------------------------
axs[0, 0].plot(t, c1, label="c1 raw")
axs[0, 0].plot(t, c1_fit, "--", label=f"c1 fit ({fit_summary['c1']['best_model']})")
axs[0, 0].set_ylabel("c1")
axs[0, 0].set_title("C1: Raw vs. Best Fit") # Individual Title
axs[0, 0].grid()
axs[0, 0].legend()

axs[0, 1].plot(t, c2, label="c2 raw")
axs[0, 1].plot(t, c2_fit, "--", label=f"c2 fit ({fit_summary['c2']['best_model']})")
axs[0, 1].set_ylabel("c2")
axs[0, 1].set_title("C2: Raw vs. Best Fit") # Individual Title
axs[0, 1].grid()
axs[0, 1].legend()

# ------------------------------------------------------------
# Row 2: c3 and c4 parameter fits
# ------------------------------------------------------------
axs[1, 0].plot(t, c3, label="c3 raw")
axs[1, 0].plot(t, c3_fit, "--", label=f"c3 fit ({fit_summary['c3']['best_model']})")
axs[1, 0].set_ylabel("c3")
axs[1, 0].set_title("C3: Raw vs. Best Fit") # Individual Title
axs[1, 0].grid()
axs[1, 0].legend()

axs[1, 1].plot(t, c4, label="c4 raw")
axs[1, 1].plot(t, c4_fit, "--", label=f"c4 fit ({fit_summary['c4']['best_model']})")
axs[1, 1].set_ylabel("c4")
axs[1, 1].set_title("C4: Raw vs. Best Fit") # Individual Title
axs[1, 1].grid()
axs[1, 1].legend()

# ------------------------------------------------------------
# Row 3: d and b parameter fits
# ------------------------------------------------------------
axs[2, 0].plot(t, d, label="d raw")
axs[2, 0].plot(t, d_fit, "--", label=f"d fit ({fit_summary['d']['best_model']})")
axs[2, 0].set_ylabel("d")
axs[2, 0].set_xlabel("Time [s]")
axs[2, 0].set_title("D: Raw vs. Best Fit") # Individual Title
axs[2, 0].grid()
axs[2, 0].legend()

axs[2, 1].plot(t, b, label="b raw")
axs[2, 1].plot(t, b_fit, "--", label=f"b fit ({fit_summary['b']['best_model']})")
axs[2, 1].set_ylabel("b")
axs[2, 1].set_xlabel("Time [s]")
axs[2, 1].set_title("B: Raw vs. Best Fit") # Individual Title
axs[2, 1].grid()
axs[2, 1].legend()

# Overall figure title for the force-model coefficient fitting results
fig.suptitle(f"Model Fitting Analysis for Components C1-C4, D, and B : Fuel level {fuel_level}%", fontsize=16)
plt.tight_layout()

# Adjust layout to make room for the overall figure title
plt.subplots_adjust(top=0.92) 



''''============================================================'''
''''---------------- Plot Force Reconstruction Results ----------'''
''''============================================================'''

# Compare measured force, raw-parameter reconstruction, and analytic-fit reconstruction
plt.figure(figsize=(12, 5))
plt.plot(t, Fx_meas, label="Fx measured", linewidth=1.5)
plt.plot(t, Fx_raw, "--", label=f"Fx reconstructed raw params  RMSE={rmse_Fx_raw:.3e}")
plt.plot(t, Fx_fitfunc, "-.", label=f"Fx reconstructed analytic fits  RMSE={rmse_Fx_fitfunc:.3e}")
plt.title(f"Force reconstruction using estimated states: Fuel level {fuel_level}%")
plt.xlabel("Time [s]")
plt.ylabel("Fx")
plt.grid()
plt.legend()

# Plot the residual between the measured force and analytic-fit reconstruction
plt.figure(figsize=(12, 5))
plt.plot(t, Fx_meas - Fx_fitfunc, label="Fx_meas - Fx_fitfunc")
plt.title(f"Force reconstruction residual using analytic parameter fits: Fuel level {fuel_level}%")
plt.xlabel("Time [s]")
plt.ylabel("Residual")
plt.grid()
plt.legend()

plt.show()
