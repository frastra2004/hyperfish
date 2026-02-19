"""
pillbox_superfish_like.py

Finite-difference axisymmetric Helmholtz eigen-solver for a pillbox cavity
(approximation in the style of Superfish). Solves for H_phi(r,z):
    (d^2/dr^2 + (1/r) d/dr + d^2/dz^2) H_phi + k^2 H_phi = 0
with Dirichlet H_phi = 0 on conducting walls. Computes eigenfrequencies,
fields, accelerating voltage on axis for beta~1, stored energy U, and R/Q.

Notes / assumptions:
- Vacuum inside cavity (eps0, mu0).
- Perfect electric conductor on all cavity walls -> Dirichlet on H_phi
  (appropriate for TM-modes formulation solved in H_phi).
- Axis r=0 handled by regularity (special FD formula).
- Uses second-order central differences.
- Outputs modes as k^2 eigenvalues; frequency f = c * sqrt(k^2) / (2*pi).
- This is a pedagogical solver — for production you should refine mesh
  and consider higher-order discretizations, geometry curvature treatment, etc.
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt 
# Physical constants
eps0 = 8.8541878128e-12
mu0  = 4*np.pi*1e-7
c0   = 1/np.sqrt(eps0*mu0)

# ---------------------------------------------------------------------
# User parameters (change these)
R = 0.05         # cavity radius [m]
L = 0.04         # cavity length [m]
Nr = 120         # number of radial grid points (including r=0)
Nz = 120         # number of axial grid points (including ends)
nev = 6          # number of eigenvalues/modes to compute
beta = 1.0       # particle beta for Vacc phase factor (use 1 for relativistic)
# ---------------------------------------------------------------------

# Build grid (cell-centered or node-centered? We'll use node-centered)
r = np.linspace(0.0, R, Nr)            # includes r=0 and r=R
z = np.linspace(0.0, L, Nz)            # includes z=0 and z=L
dr = r[1] - r[0]
dz = z[1] - z[0]
Nrn = Nr
Nzn = Nz
N = Nrn * Nzn  # total unknowns in 2D

# Index mapping: (ir, iz) -> idx = iz*Nr + ir
def idx(ir, iz):
    return iz * Nrn + ir

# Build sparse matrix A such that A * H = k^2 * H
# We'll form A = -Laplace_perp operator (so eigenproblem is A H = k^2 H)
rows = []
cols = []
data = []

# Helper to add entry
def add_entry(i, j, val):
    rows.append(i); cols.append(j); data.append(val)

# Precompute radial positions for use in 1/r term (avoid divide-by-zero)
# For r=0, we will use the regularized finite difference formula (see below)
for iz in range(Nzn):
    for ir in range(Nrn):
        i = idx(ir, iz)
        ri = r[ir]

        # boundary: Dirichlet H=0 on metal walls (r=R or z=0 or z=L)
        if (ir == Nrn-1) or (iz == 0) or (iz == Nzn-1):
            # Dirichlet: enforce H[i] = 0 -> large diag penalty (or simply set eqn as H=0)
            # Simpler: set row to identity so eigenvector value is exactly zero at boundary
            add_entry(i, i, 1.0)
            continue

        # Interior or axis points: discretize Laplacian in cylindrical coords for m=0:
        # Lap = d2/dr2 + (1/r) d/dr + d2/dz2
        # central differences:
        # d2/dr2 -> (H_{i+1} - 2 H_i + H_{i-1}) / dr^2
        # (1/r) d/dr -> (1/r) * (H_{i+1} - H_{i-1})/(2 dr)
        # d2/dz2 -> (H_{j+1} - 2 H_j + H_{j-1}) / dz^2

        # Special treatment for axis node (r=0)
        if ir == 0:
            # Use limiting form as r->0:
            # (1/r) d/dr (r dH/dr) -> 2 * d2H/dr2 at r=0 for m=0 regularity.
            # So radial part approx: 2 * (H_{ir+1} - 2 H_ir + H_{ir+2?}) / dr^2
            # Simpler and stable approach: use one-sided second-order formula:
            # H'(0) = 0 (symmetry), so impose ghost point H_{-1} = H_{1}
            # then d2/dr2 at 0: (H_1 - 2 H_0 + H_{-1})/dr^2 = 2*(H_1 - H_0)/dr^2
            # And (1/r) d/dr term tends to 0 in symmetry. So radial contribution becomes:
            # d2/dr2 + (1/r) d/dr  ≈ 2*(H_1 - H_0)/dr^2
            # Implement radial stencil accordingly.

            # indices
            i_r_plus = idx(ir+1, iz)
            # radial contribution:
            add_entry(i, i, -2.0 * (1.0 / dr**2))        # H_i coefficient comes from radial approx
            add_entry(i, i_r_plus,  2.0 * (1.0 / dr**2))

        else:
            # interior radial central difference including 1/r term
            i_r_plus  = idx(ir+1, iz)
            i_r_minus = idx(ir-1, iz)
            # coefficients for radial part
            cr_plus  = 1.0 / dr**2 + 1.0 / (2.0*dr*ri)
            cr_minus = 1.0 / dr**2 - 1.0 / (2.0*dr*ri)
            cr_center = -2.0 / dr**2
            add_entry(i, i_r_plus,  cr_plus)
            add_entry(i, i_r_minus, cr_minus)
            add_entry(i, i,        cr_center)

        # axial second derivative
        i_z_plus  = idx(ir, iz+1)
        i_z_minus = idx(ir, iz-1)
        add_entry(i, i_z_plus,  1.0 / dz**2)
        add_entry(i, i_z_minus, 1.0 / dz**2)
        add_entry(i, i,        -2.0 / dz**2)

# Assemble sparse matrix
A = sp.coo_matrix((data, (rows, cols)), shape=(N, N)).tocsr()

# For Dirichlet rows we put identity; these rows correspond to fixed zero values,
# and therefore in eigenvalue computation they will produce trivial eigenpairs.
# To avoid spurious zero rows we keep them as identity and restrict eigen-solve to interior DOFs.
# Build mask of free (non-Dirichlet) DOFs
dirichlet = np.zeros(N, dtype=bool)
for iz in range(Nzn):
    for ir in range(Nrn):
        ii = idx(ir, iz)
        if (ir == Nrn-1) or (iz == 0) or (iz == Nzn-1):
            dirichlet[ii] = True

free_dofs = np.where(~dirichlet)[0]
nfree = free_dofs.size

# Reduce the matrix to free DOFs: A_free * H_free = k^2 * H_free
A_free = A[free_dofs][:, free_dofs]

# Solve eigenvalue problem for smallest few eigenvalues (k^2)
# A_free is symmetric; use eigsh targeting smallest magnitude eigenvalues
which = 'SM'  # smallest magnitude
if nfree <= nev + 2:
    nev_calc = max(1, nfree-2)
else:
    nev_calc = nev

print("Matrix size (full):", N, "free DOFs:", nfree, "computing nev =", nev_calc)

# Use shift-invert if needed (for better accuracy near target), otherwise rely on SM
try:
    # For symmetric positive-definite-like operator, eigsh works
    eigvals, eigvecs = spla.eigsh(A_free, k=nev_calc, which=which)
except Exception as e:
    # fallback to sparse.linalg.eigs if eigsh fails
    eigvals, eigvecs = spla.eigs(A_free, k=nev_calc, which=which)
    eigvals = eigvals.real
    eigvecs = eigvecs.real

# Sort eigenpairs
perm = np.argsort(eigvals)
eigvals = eigvals[perm]
eigvecs = eigvecs[:, perm]

# Map eigenvectors back to full grid (fill Dirichlet with zeros)
modes = []
for m in range(eigvecs.shape[1]):
    Hfull = np.zeros(N, dtype=np.float64)
    Hfull[free_dofs] = eigvecs[:, m]
    modes.append(Hfull.reshape((Nzn, Nrn)))   # shape (z, r) for easier indexing

# Convert eigenvalues k^2 -> frequency (Hz), store k and f
k2_vals = eigvals
k_vals = np.sqrt(np.abs(k2_vals))
freqs = (c0 * k_vals) / (2*np.pi)

print("Eigenfrequencies (Hz) of computed modes:")
for i, f in enumerate(freqs):
    print(f"  mode {i+1}: f = {f:.6e} Hz (k^2 = {k2_vals[i]:.6e})")

# ---------------------------------------------------------------------
# Post-processing: identify TM010-like mode and compute Vacc, U, R/Q
# We'll pick the mode with maximum |E_z| on axis (E_z computed from H_phi)
# For TM (m=0), relations (frequency-domain):
#   E_r = - (i / (omega eps0)) * dH_phi/dz
#   E_z = + (i / (omega eps0)) * (1/r) * d(r H_phi)/dr
# Since eigenvectors are real-valued (up to scale), we'll compute the spatial
# patterns and compute Vacc magnitude assuming phase factor exp(i omega t).
# We'll normalize modes by stored energy U = 1 J for convenience when computing R/Q.
# ---------------------------------------------------------------------

def compute_fields_and_R_over_Q(H_grid, k):
    """
    H_grid: shape (Nz, Nr) array of H_phi on nodes
    k: wave-number (rad/m)
    returns dict with: Ez_on_axis(z), Vacc, U, R_over_Q
    """
    omega = c0 * k
    # compute E fields (real spatial patterns). We will compute E_z on axis:
    # discrete derivatives in r and z (centered)
    H = H_grid.copy()   # shape (Nz, Nr) with indices (iz, ir)

    # derivatives w.r.t r: d(r H)/dr needed for E_z
    # we compute r * H at radial nodes
    r_nodes = r.copy()
    rH = np.zeros_like(H)
    for iz in range(Nzn):
        rH[iz, :] = r_nodes * H[iz, :]

    # radial derivative d(rH)/dr using central differences (at nodes)
    dr_rH = np.zeros_like(H)
    # interior radial points
    for ir in range(1, Nrn-1):
        dr_rH[:, ir] = (rH[:, ir+1] - rH[:, ir-1]) / (2*dr)
    # axis special (one-sided): dr_rH[:,0] ≈ (rH[:,1] - rH[:,0])/dr  (consistent with earlier axis treatment)
    dr_rH[:, 0] = (rH[:, 1] - rH[:, 0]) / dr
    # at r=R boundary we don't need derivative (H=0 at boundary), but compute one-sided for completeness
    dr_rH[:, -1] = (rH[:, -1] - rH[:, -2]) / dr

    # compute E_z = (i/(omega*eps0)) * (1/r) * d(r H)/dr
    # on axis r=0 this formula must be taken with limit; E_z(0) = (i/(omega eps0)) * (d2H/dr2 at 0 times something)
    # but numerically we'll compute E_z on axis using the limit: (1/r) d(rH)/dr -> (d^2 H / dr^2 * r?) -- to avoid complications,
    # use the value at the first non-axis radial node extrapolated to axis. Simpler practical approach:
    Ez = np.zeros_like(H, dtype=np.complex128)
    for iz in range(Nzn):
        for ir in range(Nrn):
            ri = r_nodes[ir]
            if ri == 0.0:
                # use limit: (1/r) d(rH)/dr at r->0 equals 2 * dH/dr at r->0 ? but simpler: use derivative estimate from dr_rH/ r ~ dr_rH[:,1]/(r[1])
                # approximate (1/r) d(rH)/dr at axis by dr_rH[:,1] / r[1]
                val = dr_rH[iz, 1] / r_nodes[1]
            else:
                val = dr_rH[iz, ir] / ri
            Ez[iz, ir] = (1j / (omega * eps0)) * val

    # compute stored energy U = (1/4) ∫ (eps |E|^2 + mu |H|^2) dV  (time-average for phasors)
    # For axisymmetric volume element: dV = 2*pi * r * dr * dz
    Er2 = np.abs(Ez)**2  # we only used Ez; for exact U we'd include Er and H_r,H_z; but dominant Ez contribution ok for TM010
    # compute H^2 (|H_phi|^2)
    H2 = np.abs(H)**2
    # integrate
    U_e = 0.25 * eps0 * np.sum(Er2 * (2*np.pi * r_nodes) * dr * dz)
    U_m = 0.25 * mu0  * np.sum(H2  * (2*np.pi * r_nodes) * dr * dz)
    U = U_e + U_m

    # accelerating voltage for particle along center axis:
    # V_acc = ∫ E_z(0,z) * exp(i omega z / (beta c)) dz
    # take E_z at the axis index ir=0
    Ez_axis = Ez[:, 0]   # length Nz
    phase = np.exp(1j * omega * z / (beta * c0))
    Vacc = np.trapz(Ez_axis * phase, z)

    # R/Q = |V_acc|^2 / (omega * U)  (units: Ohm)
    R_over_Q = (np.abs(Vacc)**2) / (omega * U)

    return {
        'Ez_axis': Ez_axis,
        'Vacc': Vacc,
        'U': U,
        'R_over_Q': R_over_Q,
        'Ez': Ez,
        'H': H
    }

# Evaluate modes and print R/Q for each
results = []
for m_idx, H_flat in enumerate(modes):
    H_grid = H_flat.T  # currently modes stored as (z, r) via reshape above; ensure indexing consistent
    # note: earlier we used modes.append(Hfull.reshape((Nzn, Nrn))) -> shape (Nz, Nr) with [iz, ir]
    H_grid = H_flat.reshape((Nzn, Nrn))  # (iz, ir)
    k = k_vals[m_idx]
    res = compute_fields_and_R_over_Q(H_grid, k)
    # Scale the mode to impose Vacc = 3kV
    scaling_factor = 3000 / np.abs(res['Vacc'])
    H_grid *= scaling_factor
    res = compute_fields_and_R_over_Q(H_grid, k)
    results.append(res)
    f = freqs[m_idx]
    print(f"\nMode {m_idx+1}: f = {f:.6e} Hz")
    print(f"  |V_acc| = {np.abs(res['Vacc']):.6e} V")
    print(f"  U = {res['U']:.6e} J")
    print(f"  R/Q = {res['R_over_Q']:.6e} Ohm")

# Optional: plot Ez on axis for each mode
plt.figure(figsize=(8,6))
for m_idx, res in enumerate(results):       
    plt.plot(z, np.real(res['Ez_axis']), label=f'Mode {m_idx+1} (f={freqs[m_idx]:.2e} Hz)')
plt.xlabel('z (m)')
plt.ylabel('E_z on axis (V/m)')
plt.title('E_z on axis for computed modes')
plt.legend()
plt.grid()
plt.show()
# Note on normalization:
print("\nNote: eigenmodes have arbitrary scale. R/Q computed above is independent of amplitude if U and V_acc are taken from the same eigenmode.")
print("If you want R/Q in physical units independent of normalization, you can rescale mode so that U=1 J and then recompute V_acc and R/Q accordingly.")
#np.savez('pillbox_modes.npz', r=r, z=z, modes=modes, freqs=freqs, results=results)
# End of script
