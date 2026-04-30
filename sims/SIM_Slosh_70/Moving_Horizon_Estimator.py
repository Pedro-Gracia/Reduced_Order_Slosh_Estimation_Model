import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
import os

''''============================================================'''
''''--------------- CFD Data Loading Setup ----------------------'''
''''============================================================'''

# Define the fuel level associated with this estimator run
# This value is also used later when saving the estimator output file
fuel_level = 70

# Build the path to the merged CFD data file
# The data file is assumed to be located in the same folder as this script
base_dir = os.path.dirname(__file__)
com_filename = os.path.join(base_dir, "CFD_data.dat")

# Load the merged CFD dataset
# Expected columns:
#   column 0 -> simulation time
#   column 1 -> wall reaction force Fx
#   column 2 -> liquid center-of-mass displacement q
data = np.loadtxt(com_filename, comments="#", skiprows=1)

# Extract time, force measurement, and COM measurement from the merged dataset
t = data[:, 0]
Fx_meas = data[:, 1]
q_meas = data[:, 2]

# Print basic diagnostics to confirm that the CFD data was loaded correctly
print("Loaded merged dataset:")
print(f"Samples: {len(t)}")
print(f"q range [m]: {q_meas.min():.6e} to {q_meas.max():.6e}")


''''============================================================'''
''''--------------- Input Acceleration Function -----------------'''
''''============================================================'''

def a_smooth(ti):
    """
    Smooth lateral acceleration input used to drive the tank.

    This function matches the cosine-ramp acceleration profile used in the
    OpenFOAM table. The acceleration is zero before the maneuver, follows a
    smooth ramp during the forced-excitation interval, and returns to zero
    afterward.
    """
    if ti < 0.5:
        return 0.0
    elif ti <= 3.5:
        return 0.75 * 0.5 * (1 - np.cos(np.pi * (ti - 0.5) / (2.0 - 0.5)))
    else:
        return 0.0

# Evaluate the known lateral acceleration input at each CFD measurement time
a_x = np.array([a_smooth(ti) for ti in t])


''''============================================================'''
''''--------------- 2-Mode Model Dimensions ---------------------'''
''''============================================================'''

# Physical state vector:
#   x = [q1, q1d, q2, q2d]
#
# Parameter vector:
#   theta = [w1, z1, g1, w2, z2, g2, c1, c2, c3, c4, d, b]
#
# Measurement vector:
#   y = [q_total, Fx]
#
# Prior vector:
#   prior = [x0; theta]

nx = 4
nth = 12
ny = 2
nprior = nx + nth


''''============================================================'''
''''--------------- Discrete Dynamics Model ---------------------'''
''''============================================================'''

def f_disc(x, a, dt, theta):
    """
    Propagate the two-mode slosh model forward by one time step.

    Inputs:
        x     -> current modal state [q1, q1d, q2, q2d]
        a     -> known lateral acceleration input at the current time
        dt    -> time step between CFD samples
        theta -> current parameter vector

    Output:
        x_next -> predicted modal state at the next time sample
    """

    # Unpack modal states
    q1, q1d, q2, q2d = x

    # Unpack model parameters
    w1, z1, g1, w2, z2, g2, c1, c2, c3, c4, d, b = theta

    # Modal accelerations from damped second-order oscillator equations
    q1dd = -2.0 * z1 * w1 * q1d - (w1 ** 2) * q1 + g1 * a
    q2dd = -2.0 * z2 * w2 * q2d - (w2 ** 2) * q2 + g2 * a

    # Forward-Euler propagation for each modal state
    q1_next  = q1  + dt * q1d
    q1d_next = q1d + dt * q1dd
    q2_next  = q2  + dt * q2d
    q2d_next = q2d + dt * q2dd

    return np.array([q1_next, q1d_next, q2_next, q2d_next])


''''============================================================'''
''''--------------- Measurement Model ---------------------------'''
''''============================================================'''

def h_meas(x, a, theta):
    """
    Predict the measured COM displacement and wall force.

    Inputs:
        x     -> modal state [q1, q1d, q2, q2d]
        a     -> known lateral acceleration input
        theta -> current parameter vector

    Output:
        y_pred -> predicted measurement [q_total, Fx_pred]
    """

    # Unpack modal states
    q1, q1d, q2, q2d = x

    # Unpack model parameters
    w1, z1, g1, w2, z2, g2, c1, c2, c3, c4, d, b = theta

    # COM displacement is modeled as the sum of the two modal coordinates
    q_total = q1 + q2

    # Wall force model using modal states, acceleration feedthrough, and bias
    Fx_pred = c1*q1 + c2*q1d + c3*q2 + c4*q2d + d*a + b

    return np.array([q_total, Fx_pred])


''''============================================================'''
''''--------------- Augmented Dynamics Jacobian -----------------'''
''''============================================================'''

def A_aug(x, a, dt, theta):
    """
    Local augmented dynamics Jacobian.

    This matrix linearizes the augmented propagation model:
        X_{k+1} = F(X_k, a_k)

    where:
        X = [x; theta] in R^16

    The physical states propagate through the discrete dynamics, while the
    parameters are modeled as locally constant over one propagation step.
    """

    # Unpack modal states
    q1, q1d, q2, q2d = x

    # Unpack model parameters
    w1, z1, g1, w2, z2, g2, c1, c2, c3, c4, d, b = theta

    # Initialize augmented Jacobian as identity
    A = np.eye(nx + nth)

    # ------------------------------------------------------------
    # State-to-state block
    # ------------------------------------------------------------
    A[0, 1] = dt
    A[1, 0] = -dt * (w1 ** 2)
    A[1, 1] = 1.0 - dt * (2.0 * z1 * w1)

    A[2, 3] = dt
    A[3, 2] = -dt * (w2 ** 2)
    A[3, 3] = 1.0 - dt * (2.0 * z2 * w2)

    # ------------------------------------------------------------
    # State-to-parameter block
    # ------------------------------------------------------------

    # q1_next has no direct theta dependence in this Euler update

    # q1d_next parameter sensitivities
    A[1, 4] = dt * (-2.0 * z1 * q1d - 2.0 * w1 * q1)  # d/dw1
    A[1, 5] = dt * (-2.0 * w1 * q1d)                  # d/dz1
    A[1, 6] = dt * a                                  # d/dg1

    # q2d_next parameter sensitivities
    A[3, 7] = dt * (-2.0 * z2 * q2d - 2.0 * w2 * q2) # d/dw2
    A[3, 8] = dt * (-2.0 * w2 * q2d)                 # d/dz2
    A[3, 9] = dt * a                                 # d/dg2

    # Theta-to-theta block remains identity because parameters are held
    # constant over one local propagation step
    return A


''''============================================================'''
''''--------------- Augmented Measurement Jacobian --------------'''
''''============================================================'''

def H_aug(x, a, theta):
    """
    Local augmented measurement Jacobian.

    This matrix linearizes the measurement model:
        y = h(X, a)

    where:
        X = [x; theta] in R^16

    The first measurement is COM displacement and the second measurement is
    wall reaction force.
    """

    # Unpack modal states
    q1, q1d, q2, q2d = x

    # Unpack model parameters
    w1, z1, g1, w2, z2, g2, c1, c2, c3, c4, d, b = theta

    # Initialize measurement Jacobian
    H = np.zeros((ny, nx + nth))

    # ------------------------------------------------------------
    # COM measurement:
    #   q = q1 + q2
    # ------------------------------------------------------------
    H[0, 0] = 1.0
    H[0, 2] = 1.0

    # ------------------------------------------------------------
    # Force measurement:
    #   Fx = c1*q1 + c2*q1d + c3*q2 + c4*q2d + d*a + b
    # ------------------------------------------------------------
    H[1, 0] = c1
    H[1, 1] = c2
    H[1, 2] = c3
    H[1, 3] = c4

    H[1, 10] = q1
    H[1, 11] = q1d
    H[1, 12] = q2
    H[1, 13] = q2d
    H[1, 14] = a
    H[1, 15] = 1.0

    return H


''''============================================================'''
''''--------------- MHE Window and Update Settings --------------'''
''''============================================================'''

# Number of samples in each moving horizon window
Nh = 50

# Number of samples advanced between consecutive MHE solves
N_update = 5


''''============================================================'''
''''--------------- Noise and Covariance Tuning -----------------'''
''''============================================================'''

# Measurement covariance matrix
# The first entry weights the COM displacement residual.
# The second entry weights the wall force residual.
R_meas = np.diag([
    (5e-4)**2,   # variance of q residual
    (1.0)**2     # variance of Fx residual
])

# Dynamic-defect covariance for the physical states
# These values control how strongly the optimizer enforces the reduced-order
# model dynamics between consecutive nodes.
Q_dyn = np.diag([
    (2e-4)**2,
    (2e-3)**2,
    (2e-4)**2,
    (2e-3)**2
])

# Parameter-drift inflation between windows
# This prevents the arrival covariance from becoming too confident and allows
# the estimated parameters to slowly adapt from one MHE window to the next.
Q_theta_drift = np.diag([
    2e-4, 1e-5, 5e-5,
    2e-4, 1e-5, 5e-5,
    1e-2, 1e-2, 1e-2, 1e-2,
    1e-4, 1e-4
])


''''============================================================'''
''''--------------- Initial Parameter Guess and Bounds ----------'''
''''============================================================'''

# Initial parameter guess used to initialize the first MHE window
theta_guess = np.array([
    2.6,    # w1
    0.08,   # z1
    0.18,   # g1
    2.9,    # w2
    0.18,   # z2
    0.09,   # g2
    290.0,  # c1
    290.0,  # c2
    -170.0, # c3
    -260.0, # c4
    2.5,    # d
    0.3     # b
], dtype=float)

# Lower bounds for the estimated parameter vector
theta_lb = np.array([
    0.1, 1e-4, 0.0,
    0.1, 1e-4, 0.0,
   -500.0, -500.0, -500.0, -500.0,
   -20.0, -20.0
])

# Upper bounds for the estimated parameter vector
theta_ub = np.array([
    20.0, 1.0, 10.0,
    20.0, 1.0, 10.0,
    500.0, 500.0, 500.0, 500.0,
    20.0, 20.0
])


''''============================================================'''
''''--------------- Decision Vector Utilities -------------------'''
''''============================================================'''

def pack_decision(x_seq, theta):
    """
    Pack the MHE decision variables into one optimization vector.

    Decision vector structure:
        z = [x_0, x_1, ..., x_N, theta]
    """
    return np.concatenate([x_seq.reshape(-1), theta])


def unpack_decision(z, N_nodes):
    """
    Unpack the MHE decision vector into state history and parameter vector.

    Inputs:
        z       -> optimizer decision vector
        N_nodes -> number of nodes in the current MHE window

    Outputs:
        x_seq -> estimated state sequence over the window
        theta -> estimated parameter vector over the window
    """
    x_seq = z[:N_nodes*nx].reshape((N_nodes, nx))
    theta = z[N_nodes*nx:]
    return x_seq, theta


def safe_chol_inv(P):
    """
    Compute a numerically safe Cholesky factor of P^{-1}.

    This is used to convert covariance-weighted residuals into whitened
    residuals for least_squares.
    """

    # Small regularization to keep the covariance symmetric and invertible
    eps = 1e-10
    Preg = 0.5 * (P + P.T) + eps * np.eye(P.shape[0])

    # Use pseudo-inverse for numerical robustness
    Pinv = np.linalg.pinv(Preg)

    # Return lower Cholesky factor of the inverse covariance
    return np.linalg.cholesky(Pinv + 1e-12*np.eye(P.shape[0]))


def build_initial_window_guess(t_win, q_win, qdot_win, mu_prior):
    """
    Build the initial guess for the current MHE window.

    Inputs:
        t_win     -> time samples in the current window
        q_win     -> measured COM displacement in the current window
        qdot_win  -> numerical derivative of measured COM displacement
        mu_prior  -> current arrival prior mean [x0_prior; theta_prior]

    Outputs:
        x_guess          -> initial state trajectory guess
        theta_prior.copy -> initial parameter guess for this window
    """

    # Split the arrival mean into initial physical state and parameter guess
    x0_prior = mu_prior[:nx]
    theta_prior = mu_prior[nx:]

    # Initialize the state trajectory guess
    x_guess = np.zeros((len(t_win), nx))
    x_guess[0] = x0_prior.copy()

    # Roll out the reduced-order model from the arrival-prior state
    for i in range(len(t_win) - 1):
        dt = t_win[i+1] - t_win[i]
        a = a_win_global_lookup(t_win[i])
        x_guess[i+1] = f_disc(x_guess[i], a, dt, theta_prior)

    # Blend the model rollout with measured COM information
    # This gives the optimizer a better starting point without fully forcing
    # the modal states to equal the measured total COM motion.
    x_guess[:, 0] = 0.7*x_guess[:, 0] + 0.3*(0.7*q_win)
    x_guess[:, 2] = 0.7*x_guess[:, 2] + 0.3*(0.3*q_win)
    x_guess[:, 1] = 0.7*x_guess[:, 1] + 0.3*(0.7*qdot_win)
    x_guess[:, 3] = 0.7*x_guess[:, 3] + 0.3*(0.3*qdot_win)

    return x_guess, theta_prior.copy()


def a_win_global_lookup(ti):
    """
    Interpolate the known acceleration input at an arbitrary time value.
    """
    return np.interp(ti, t, a_x)


''''============================================================'''
''''--------------- MHE Residual Function -----------------------'''
''''============================================================'''

def mhe_residuals(z, t_win, q_win, Fx_win, a_win, mu_prior, P_prior):
    """
    Assemble the full MHE residual vector for least_squares.

    Decision variables:
        z = [x_0, x_1, ..., x_N, theta]

    Residual blocks:
        1) arrival prior residual
        2) dynamic-defect residuals
        3) measurement residuals

    The arrival prior is joint on:
        [x_0; theta]
    """

    # Number of nodes in the current MHE window
    N_nodes = len(t_win)

    # Recover state sequence and parameter vector from optimizer variables
    x_seq, theta = unpack_decision(z, N_nodes)

    # Initialize residual list
    res = []

    # ------------------------------------------------------------
    # Arrival residual
    # ------------------------------------------------------------
    # This penalizes deviation from the propagated prior mean/covariance.
    z_prior_vec = np.concatenate([x_seq[0], theta])
    L_prior = safe_chol_inv(P_prior)
    res.extend(L_prior @ (z_prior_vec - mu_prior))

    # ------------------------------------------------------------
    # Dynamics residuals
    # ------------------------------------------------------------
    # These penalize deviations between consecutive MHE states and the
    # reduced-order model prediction.
    LQ = safe_chol_inv(Q_dyn)

    for i in range(N_nodes - 1):
        dt = t_win[i+1] - t_win[i]
        x_pred = f_disc(x_seq[i], a_win[i], dt, theta)
        defect = x_seq[i+1] - x_pred
        res.extend(LQ @ defect)

    # ------------------------------------------------------------
    # Measurement residuals
    # ------------------------------------------------------------
    # These compare the predicted measurements against the CFD-derived
    # COM displacement and wall force data.
    LR = safe_chol_inv(R_meas)

    for i in range(N_nodes):
        y_pred = h_meas(x_seq[i], a_win[i], theta)
        y_meas_i = np.array([q_win[i], Fx_win[i]])
        res.extend(LR @ (y_pred - y_meas_i))

    return np.array(res)


''''============================================================'''
''''--------------- Initial Arrival Mean and Covariance ---------'''
''''============================================================'''

# Estimate COM velocity from the measured COM displacement
qdot_guess = np.gradient(q_meas, t)

# Initial arrival mean:
# The measured COM is split between the two modes using a 70/30 split.
mu_prior = np.array([
    0.7*q_meas[0],
    0.7*qdot_guess[0],
    0.3*q_meas[0],
    0.3*qdot_guess[0],
    *theta_guess
], dtype=float)

# Initial arrival covariance:
# The first four entries correspond to the initial modal states.
# The remaining entries correspond to the initial parameter uncertainty.
P_prior = np.diag([
    1e-5, 1e-3, 1e-5, 1e-3,   # x0 prior covariance
    0.2**2, 0.03**2, 0.05**2, # w1,z1,g1
    0.2**2, 0.03**2, 0.05**2, # w2,z2,g2
    10.0**2, 10.0**2, 10.0**2, 10.0**2,
    0.5**2, 0.5**2
])


''''============================================================'''
''''--------------- Estimator Storage Arrays --------------------'''
''''============================================================'''

# Number of CFD measurement samples
N = len(t)

# Estimated modal states at each time sample
x_hat_hist = np.full((N, nx), np.nan)

# Estimated parameter vector at each time sample
theta_hist = np.full((N, nth), np.nan)

# Predicted measurement history [q_est, Fx_est]
y_hat_hist = np.full((N, ny), np.nan)

# Posterior covariance history for diagnostics
P_post_hist = np.full((N, nx+nth, nx+nth), np.nan)


''''============================================================'''
''''--------------- Moving Horizon Estimation Loop --------------'''
''''============================================================'''

for k_end in range(Nh, N, N_update):

    # Define the start and end index of the current MHE window
    k_start = k_end - Nh

    # Extract current MHE window data
    t_win = t[k_start:k_end+1]
    q_win = q_meas[k_start:k_end+1]
    Fx_win = Fx_meas[k_start:k_end+1]
    a_win = a_x[k_start:k_end+1]
    qdot_win = qdot_guess[k_start:k_end+1]

    # Build initial guess for the current MHE solve
    x_guess, theta_guess_local = build_initial_window_guess(t_win, q_win, qdot_win, mu_prior)
    z0 = pack_decision(x_guess, theta_guess_local)

    # Number of nodes in current MHE window
    N_nodes = len(t_win)

    # Bounds for modal state sequence over the window
    x_lb = np.tile(np.array([-0.05, -1.0, -0.05, -1.0]), N_nodes)
    x_ub = np.tile(np.array([ 0.05,  1.0,  0.05,  1.0]), N_nodes)

    # Full lower and upper bounds for optimizer decision vector
    z_lb = np.concatenate([x_lb, theta_lb])
    z_ub = np.concatenate([x_ub, theta_ub])

    # Solve the nonlinear least-squares MHE problem
    result = least_squares(
        mhe_residuals,
        z0,
        bounds=(z_lb, z_ub),
        args=(t_win, q_win, Fx_win, a_win, mu_prior, P_prior),
        verbose=0,
        max_nfev=40
    )

    # Extract optimized state sequence and parameter vector
    z_opt = result.x
    x_opt, theta_opt = unpack_decision(z_opt, N_nodes)

    # Store predicted measurements, modal states, and parameters over this window
    for i in range(N_nodes):
        y_hat_hist[k_start+i] = h_meas(x_opt[i], a_win[i], theta_opt)
        x_hat_hist[k_start+i] = x_opt[i]
        theta_hist[k_start+i] = theta_opt

    # ------------------------------------------------------------
    # Posterior covariance approximation from J^T J
    # ------------------------------------------------------------
    # The least-squares Jacobian provides a local approximation of the
    # inverse posterior covariance.
    J = result.jac
    JTJ = J.T @ J

    # Regularized inverse for numerical stability
    Sigma_post = np.linalg.pinv(JTJ + 1e-10*np.eye(JTJ.shape[0]))

    # Select the state block that will become the arrival state for the
    # start of the next MHE window
    shift = min(N_update, N_nodes - 1)

    idx_x_shift = np.arange(shift*nx, (shift+1)*nx)
    idx_theta = np.arange(N_nodes*nx, N_nodes*nx + nth)
    idx_joint = np.concatenate([idx_x_shift, idx_theta])

    # Extract joint posterior covariance and mean for [x_shift; theta]
    P_joint_post = Sigma_post[np.ix_(idx_joint, idx_joint)]
    mu_joint_post = np.concatenate([x_opt[shift], theta_opt])

    # ------------------------------------------------------------
    # Parameter-drift covariance inflation
    # ------------------------------------------------------------
    # Symmetrize covariance and inflate the parameter block so the next
    # MHE window can continue adapting the time-varying parameters.
    P_joint_post = 0.5*(P_joint_post + P_joint_post.T)
    P_joint_post[nx:, nx:] += N_update * Q_theta_drift

    # Update arrival prior for the next MHE window
    mu_prior = mu_joint_post.copy()
    P_prior = P_joint_post.copy()

    # Store posterior covariance at the shifted endpoint for diagnostics
    end_idx = k_start + shift
    P_post_hist[end_idx] = P_joint_post

    # Print progress occasionally during the MHE solve loop
    if (k_end // N_update) % 20 == 0:
        print(f"MHE with covariance update solved window ending at t = {t[k_end]:.2f} s")


''''============================================================'''
''''--------------- Fill Leading NaN Values ---------------------'''
''''============================================================'''

# Before the first solved MHE window, the storage arrays contain NaNs.
# Fill these leading values using the first valid estimator result so the
# plotting and saved arrays remain continuous.
first_valid = np.where(~np.isnan(x_hat_hist[:, 0]))[0]
if len(first_valid) > 0:
    i0 = first_valid[0]
    x_hat_hist[:i0] = x_hat_hist[i0]
    theta_hist[:i0] = theta_hist[i0]
    y_hat_hist[:i0] = y_hat_hist[i0]


''''============================================================'''
''''--------------- Extract Estimated States and Parameters -----'''
''''============================================================'''

# Estimated modal states
q1_est = x_hat_hist[:, 0]
q1d_est = x_hat_hist[:, 1]
q2_est = x_hat_hist[:, 2]
q2d_est = x_hat_hist[:, 3]

# Estimated modal parameters
w1_est = theta_hist[:, 0]
z1_est = theta_hist[:, 1]
g1_est = theta_hist[:, 2]
w2_est = theta_hist[:, 3]
z2_est = theta_hist[:, 4]
g2_est = theta_hist[:, 5]

# Estimated force-model parameters
c1_est = theta_hist[:, 6]
c2_est = theta_hist[:, 7]
c3_est = theta_hist[:, 8]
c4_est = theta_hist[:, 9]
d_est  = theta_hist[:, 10]
b_est  = theta_hist[:, 11]

# Estimated measurement outputs
q_est = y_hat_hist[:, 0]
Fx_est = y_hat_hist[:, 1]


''''============================================================'''
''''--------------- Save Estimator Results ----------------------'''
''''============================================================'''

# Save estimated parameter histories, modal states, COM fit, and force fit
# These outputs are used later for parameter fitting and force reconstruction.
np.savez(f"Estimator_results_{fuel_level}.npz",
         t_30=t,
         w1_30=w1_est,
         z1_30=z1_est,
         g1_30=g1_est,
         w2_30=w2_est,
         z2_30=z2_est,
         g2_30=g2_est,
         c1_30=c1_est,
         c2_30=c2_est,
         c3_30=c3_est,
         c4_30=c4_est,
         d_30 =d_est,
         b_30 =b_est,
         q_30 = q_est,
         Fx_30 = Fx_est,
         q1_30 = q1_est,
         q1d_30 = q1d_est,
         q2_30 = q2_est,  
         q2d_30 = q2d_est
)


''''============================================================'''
''''--------------- Summary Printout ----------------------------'''
''''============================================================'''

# Select the forced-excitation region using a threshold on the input acceleration
a_thresh = 0.05 * np.max(np.abs(a_x))
idx = np.abs(a_x) > a_thresh

# Print mean estimated parameters over the selected excitation interval
print("\n===== COVARIANCE-UPDATED MHE STEADY PARAMETER MEANS (t > 5 s) =====")
print(f"w1 = {np.nanmean(w1_est[idx]):.6f}")
print(f"z1 = {np.nanmean(z1_est[idx]):.6f}")
print(f"g1 = {np.nanmean(g1_est[idx]):.6f}") 
print(f"w2 = {np.nanmean(w2_est[idx]):.6f}")
print(f"z2 = {np.nanmean(z2_est[idx]):.6f}")
print(f"g2 = {np.nanmean(g2_est[idx]):.6f}")
print(f"c1 = {np.nanmean(c1_est[idx]):.6f}")
print(f"c2 = {np.nanmean(c2_est[idx]):.6f}")
print(f"c3 = {np.nanmean(c3_est[idx]):.6f}")
print(f"c4 = {np.nanmean(c4_est[idx]):.6f}")
print(f"d  = {np.nanmean(d_est[idx]):.6f}")
print(f"b  = {np.nanmean(b_est[idx]):.6f}")


''''============================================================'''
''''--------------- COM and Force Fit Plots ---------------------'''
''''============================================================'''

plt.figure(figsize=(12, 8))

# Plot measured COM displacement against MHE reconstructed COM displacement
plt.subplot(2, 1, 1)
plt.plot(t, q_meas, label="CFD q")
plt.plot(t, q_est, "--", label="MHE q (cov-updated)")
plt.axvline(5.0, color='k', linestyle='--', linewidth=1)
plt.ylabel("q [m]")
plt.title("2-Mode MHE with Prior/Posterior Covariance Update: COM Fit")
plt.grid()
plt.legend()

# Plot measured wall force against MHE reconstructed wall force
plt.subplot(2, 1, 2)
plt.plot(t, Fx_meas, label="CFD Fx")
plt.plot(t, Fx_est, "--", label="MHE Fx (cov-updated)")
plt.axvline(5.0, color='k', linestyle='--', linewidth=1)
plt.xlabel("Time [s]")
plt.ylabel("Fx [N]")
plt.title("2-Mode MHE with Prior/Posterior Covariance Update: Force Fit")
plt.grid()
plt.legend()

plt.tight_layout()
plt.show()


''''============================================================'''
''''--------------- Estimated Parameter History Plots -----------'''
''''============================================================'''

plt.figure(figsize=(12, 12))

# Mode 1 natural frequency history
plt.subplot(5, 2, 1)
plt.plot(t, w1_est)
plt.ylabel("w1")
plt.grid()

# Mode 1 damping ratio history
plt.subplot(5, 2, 2)
plt.plot(t, z1_est)
plt.ylabel("z1")
plt.grid()

# Mode 1 acceleration gain history
plt.subplot(5, 2, 3)
plt.plot(t, g1_est)
plt.ylabel("g1")
plt.grid()

# Mode 2 natural frequency history
plt.subplot(5, 2, 4)
plt.plot(t, w2_est)
plt.ylabel("w2")
plt.grid()

# Mode 2 damping ratio history
plt.subplot(5, 2, 5)
plt.plot(t, z2_est)
plt.ylabel("z2")
plt.grid()

# Mode 2 acceleration gain history
plt.subplot(5, 2, 6)
plt.plot(t, g2_est)
plt.ylabel("g2")
plt.grid()

# Force coefficients associated with mode 1
plt.subplot(5, 2, 7)
plt.plot(t, c1_est, label="c1")
plt.plot(t, c2_est, label="c2")
plt.ylabel("c1,c2")
plt.grid()
plt.legend()

# Force coefficients associated with mode 2
plt.subplot(5, 2, 8)
plt.plot(t, c3_est, label="c3")
plt.plot(t, c4_est, label="c4")
plt.ylabel("c3,c4")
plt.grid()
plt.legend()

# Acceleration feedthrough and force bias histories
plt.subplot(5, 2, 9)
plt.plot(t, d_est, label="d")
plt.plot(t, b_est, label="b")
plt.ylabel("d,b")
plt.xlabel("Time [s]")
plt.grid()
plt.legend()

plt.tight_layout()
plt.show()
