"""
pillbox_superfish_like_with_Ezero_and_target_freq.py

Finite-difference axisymmetric Helmholtz eigen-solver for a pillbox cavity
(Superfish-like). Added features:
 - pick eigenmode nearest to a user-specified target frequency (Hz)
 - scale mode so peak on-axis E_z magnitude == Ezero (V/m)
 - plot E_z magnitude along axis

Usage: edit user parameters below and run with Python 3.x.

Dependencies:
  numpy, scipy, matplotlib

Notes:
 - This is a pedagogical FD solver for TM modes (H_phi scalar).
 - The R/Q computed is independent of amplitude; Vacc scales with amplitude.
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

# ----------------------------
# Physical constants
# ----------------------------
eps0 = 8.8541878128e-12
mu0  = 4*np.pi*1e-7
c0   = 1/np.sqrt(eps0*mu0)

# ----------------------------
# User parameters
# ----------------------------
R = 0.05         # cavity radius [m]
L = 0.04         # cavity length [m]
Nr = 160         # radial grid points including r=0 and r=R
Nz = 160         # axial grid points including z=0 and z=L
nev = 8          # number of eigenvalues to compute
beta = 1.0       # particle beta for Vacc phase factor
target_freq = 3.0e9   # target frequency in Hz (set to None to pick lowest mode)
Ezero = 5e6           # desired peak on-axis Ez magnitude in V/m (set to None to leave mode unscaled)
# ----------------------------

# build grid
r = np.linspace(0.0, R, Nr)
z = np.linspace(0.0, L, Nz)
dr = r[1] - r[0]
dz = z[1] - z[0]
Nrn = Nr
Nzn = Nz
N = Nrn * Nzn

def idx(ir, iz):
    return iz * Nrn + ir

# assemble sparse matrix A such that A * H = k^2 * H
rows = []
cols = []
data = []

def add_entry(i, j, val):
    rows.append(i); cols.append(j); data.append(val)

for iz in range(Nzn):
    for ir in range(Nrn):
        i = idx(ir, iz)
        ri = r[ir]

        # Dirichlet: H=0 on metallic walls at r=R, z=0, z=L
        if (ir == Nrn-1) or (iz == 0) or (iz == Nzn-1):
            add_entry(i, i, 1.0)
            continue

        # interior or axis
        if ir == 0:
            # axis special: radial radial part approx -> 2*(H_1 - H_0)/dr^2 for d2/dr2 + (1/r) d/dr limit
            add_entry(i, i, -2.0 / dr**2)
            add_entry(i, idx(ir+1, iz), 2.0 / dr**2)
        else:
            i_r_plus  = idx(ir+1, iz)
            i_r_minus = idx(ir-1, iz)
            cr_plus  = 1.0 / dr**2 + 1.0 / (2.0*dr*ri)
            cr_minus = 1.0 / dr**2 - 1.0 / (2.0*dr*ri)
            cr_center = -2.0 / dr**2
            add_entry(i, i_r_plus,  cr_plus)
            add_entry(i, i_r_minus, cr_minus)
            add_entry(i, i,        cr_center)

        # axial second derivative
        add_entry(i, idx(ir, iz+1), 1.0 / dz**2)
        add_entry(i, idx(ir, iz-1), 1.0 / dz**2)
        add_entry(i, i,            -2.0 / dz**2)

A = sp.coo_matrix((data, (rows, cols)), shape=(N, N)).tocsr()

# identify Dirichlet rows and free DOFs
dirichlet = np.zeros(N, dtype=bool)
for iz in range(Nzn):
    for ir in range(Nrn):
        ii = idx(ir, iz)
        if (ir == Nrn-1) or (iz == 0) or (iz == Nzn-1):
            dirichlet[ii] = True
free_dofs = np.where(~dirichlet)[0]
nfree = free_dofs.size
A_free = A[free_dofs][:, free_dofs]

# eigen-solve (smallest magnitude eigenvalues)
which = 'SM'
nev_calc = min(nev, max(1, nfree-2))
print(f"Matrix full size: {N}, free DOFs: {nfree}, computing {nev_calc} eigenpairs...")

try:
    eigvals, eigvecs = spla.eigsh(A_free, k=nev_calc, which=which)
except Exception as e:
    eigvals, eigvecs = spla.eigs(A_free, k=nev_calc, which=which)
    eigvals = eigvals.real
    eigvecs = eigvecs.real

perm = np.argsort(eigvals)
eigvals = eigvals[perm]
eigvecs = eigvecs[:, perm]

k2_vals = eigvals
# remove small negative numerical noise
k2_vals = np.maximum(k2_vals, 0.0)
k_vals = np.sqrt(k2_vals)
freqs = (c0 * k_vals) / (2*np.pi)

print("Computed eigenfrequencies (Hz):")
for i, f in enumerate(freqs):
    print(f"  mode {i+1}: {f:.9e}")

# map eigenvectors back to full grid and store as (iz, ir)
modes = []
for m in range(eigvecs.shape[1]):
    Hfull = np.zeros(N, dtype=np.float64)
    Hfull[free_dofs] = eigvecs[:, m]
    modes.append(Hfull.reshape((Nzn, Nrn)))   # (iz, ir)

# helper: compute Ez, Vacc, U, R/Q for a given H_grid and k
def compute_fields_and_R_over_Q(H_grid, k, do_full_energy=True):
    omega = c0 * k
    H = H_grid.copy()   # shape (Nz, Nr)
    r_nodes = r.copy()
    # compute r*H and its radial derivative
    rH = (r_nodes[np.newaxis, :] * H)
    dr_rH = np.zeros_like(H)
    # interior radial central differences
    for ir in range(1, Nrn-1):
        dr_rH[:, ir] = (rH[:, ir+1] - rH[:, ir-1]) / (2*dr)
    # axis and outer
    dr_rH[:, 0] = (rH[:, 1] - rH[:, 0]) / dr
    dr_rH[:, -1] = (rH[:, -1] - rH[:, -2]) / dr

    # compute Ez (complex phasor): Ez = (i/(omega*eps0)) * (1/r) * d(rH)/dr
    Ez = np.zeros_like(H, dtype=np.complex128)
    for iz in range(Nzn):
        for ir in range(Nrn):
            ri = r_nodes[ir]
            if ri == 0.0:
                # approximate using first radial node derivative scaled by r[1]
                val = dr_rH[iz, 1] / r_nodes[1]
            else:
                val = dr_rH[iz, ir] / ri
            Ez[iz, ir] = (1j / (omega * eps0)) * val

    # approximate total fields for stored energy: compute only Ez and H_phi contributions
    Er2 = np.abs(Ez)**2
    H2  = np.abs(H)**2
    # axisymmetric volume element: 2*pi*r dr dz
    integrand_E = Er2 * (2*np.pi * r_nodes[np.newaxis, :])
    integrand_H = H2  * (2*np.pi * r_nodes[np.newaxis, :])
    Ue = 0.25 * eps0 * np.sum(integrand_E) * dr * dz
    Um = 0.25 * mu0  * np.sum(integrand_H) * dr * dz
    U = Ue + Um

    # accelerating voltage along axis (take Ez at ir=0)
    Ez_axis = Ez[:, 0]   # shape Nz
    phase = np.exp(1j * omega * z / (beta * c0))
    Vacc = np.trapz(Ez_axis * phase, z)

    R_over_Q = (np.abs(Vacc)**2) / (omega * U) if U != 0 else np.nan

    return {
        'Ez': Ez,
        'Ez_axis': Ez_axis,
        'Vacc': Vacc,
        'U': U,
        'R_over_Q': R_over_Q
    }

# choose mode: nearest to target_freq if provided, otherwise pick lowest mode
if target_freq is None:
    chosen_index = 0
else:
    # compute absolute difference and pick smallest
    diffs = np.abs(freqs - target_freq)
    chosen_index = int(np.argmin(diffs))
print(f"Chosen mode index: {chosen_index+1} with frequency {freqs[chosen_index]:.9e} Hz (target {target_freq})")

# compute fields for chosen mode (unscaled)
H_grid = modes[chosen_index]   # shape (Nz, Nr)
k_chosen = k_vals[chosen_index]
res = compute_fields_and_R_over_Q(H_grid, k_chosen)

print("\nBefore scaling:")
print(f"  Mode freq = {freqs[chosen_index]:.6e} Hz")
print(f"  |V_acc| (unscaled) = {np.abs(res['Vacc']):.6e} (arbitrary units)")
print(f"  U (unscaled) = {res['U']:.6e} J")
print(f"  R/Q (geometry) = {res['R_over_Q']:.6e} Ohm")

# scaling: if Ezero provided, scale fields so that peak |Ez_axis| == Ezero
scale_factor = 1.0
if Ezero is not None:
    Ez_axis_unscaled = res['Ez_axis']
    peak_Ez = np.max(np.abs(Ez_axis_unscaled))
    if peak_Ez == 0:
        raise RuntimeError("Unscaled mode has zero on-axis Ez; cannot scale to Ezero.")
    scale_factor = Ezero / peak_Ez
    print(f"\nScaling mode by factor {scale_factor:.6e} so peak |Ez_axis| == {Ezero:.6e} V/m")
    # apply scaling to H (and thereby Ez when recomputed)
    H_grid_scaled = H_grid * scale_factor
    res_scaled = compute_fields_and_R_over_Q(H_grid_scaled, k_chosen)
else:
    H_grid_scaled = H_grid.copy()
    res_scaled = res

# final printed results
omega = c0 * k_chosen
print("\nAfter scaling (if applied):")
print(f"  Mode freq = {freqs[chosen_index]:.6e} Hz")
print(f"  peak |Ez_axis| = {np.max(np.abs(res_scaled['Ez_axis'])):.6e} V/m")
print(f"  |V_acc| = {np.abs(res_scaled['Vacc']):.6e} V")
print(f"  U = {res_scaled['U']:.6e} J")
print(f"  R/Q = {res_scaled['R_over_Q']:.6e} Ohm")

# Plot Ez along axis (magnitude) vs z
Ez_axis_plot = res_scaled['Ez_axis']
plt.figure(figsize=(8,4))
plt.plot(z, np.abs(Ez_axis_plot), '-o', markersize=3, label='|E_z(on-axis)|')
plt.xlabel('z (m)')
plt.ylabel(r'$|E_z|$ (V/m)')
plt.title(f'On-axis $E_z$ magnitude — mode {chosen_index+1}, f={freqs[chosen_index]:.3e} Hz')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Optionally: plot 2D Ez map (magnitude) for visualization
# Uncomment below to plot a 2D color image of |Ez|(r,z)
Ez2D = np.abs(res_scaled['Ez'])
plt.figure(figsize=(6,5))
rr, zz = np.meshgrid(r, z)
plt.pcolormesh(rr, zz, Ez2D, shading='auto')
plt.colorbar(label='|E_z| (V/m)')
plt.xlabel('r (m)')
plt.ylabel('z (m)')
plt.title('Ez magnitude (r,z)')
plt.tight_layout()
plt.show()

# End of script
