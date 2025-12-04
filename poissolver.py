import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

eps0 = 8.854187817e-12

def assemble_poisson_cylindrical(R, L, Nr, Nz):
    """
    Build sparse matrix A for Laplacian in cylindrical coords (axisymmetric)
    on a grid r = [0..R], z = [0..L] with Nr x Nz grid points (including boundaries).
    Unknowns are interior points only (exclude Dirichlet-boundary nodes).
    Returns: A (csr), mapping info (r_grid, z_grid, idx2ij, ij2idx), hx, hz
    """
    # grid including boundaries
    r = np.linspace(0, R, Nr)
    z = np.linspace(0, L, Nz)
    dr = r[1] - r[0]
    dz = z[1] - z[0]

    # interior indices: i = 1..Nr-2 (r), j = 1..Nz-2 (z)
    nr = Nr - 2
    nz = Nz - 2
    N = nr * nz
    if N <= 0:
        raise ValueError("Need at least 3 points in each direction (including boundaries).")

    # helper maps
    def ij_to_k(i_r, j_z):  # i_r = 0..nr-1 ; j_z = 0..nz-1
        return j_z * nr + i_r

    # assemble lists for sparse matrix
    rows = []
    cols = []
    data = []

    for j in range(nz):
        for i in range(nr):
            k = ij_to_k(i, j)
            # physical coordinates of interior node
            rp = r[i+1]   # interior r index offset by 1 (since boundaries at 0 and Nr-1)
            # finite-difference coefficients for (1/r) d/dr (r dφ/dr) + d2φ/dz2
            # center coefficient
            Cc = 0.0

            # r-direction neighbors
            # left neighbor (r_{i-1})  -- careful when i==0 (adjacent to axis r=0)
            if i == 0:
                # node sits adjacent to r=0; apply symmetry ∂φ/∂r = 0 implying ghost value φ(-dr)=φ(+dr).
                # Implement by using one-sided stencil: approximate radial second derivative as
                # (φ_{i+1} - φ_i) * 2 / dr^2  (from symmetry), and the (1/r) term handled via limit r->0:
                # At r=0 we can use (1/r)d/dr(r dφ/dr) -> d^2φ/dr^2 at r=0.
                # We'll implement a compact form: use central-like: coeff_center += -2/dr^2; coeff_right += 2/dr^2
                a_center_r = -2.0 / (dr*dr)
                a_right_r  =  2.0 / (dr*dr)
                Cc += a_center_r
                # right neighbor exists (i+1)
                rows.append(k); cols.append(ij_to_k(i+1, j)); data.append(a_right_r)
            else:
                # general interior r not axis:
                a_left_r  = (1.0 / (2*dr)) * (1.0/(rp - dr/2.0)) + 1.0/(dr*dr)  # derived from discretization
                a_right_r = -(1.0 / (2*dr)) * (1.0/(rp + dr/2.0)) + 1.0/(dr*dr)
                # but simpler (more stable) approach: use standard central differences for r-derivatives:
                # (1/r) d/dr (r dφ/dr) ≈ (1/rp) * ( (rp+1/2) (φ_{i+1}-φ_i)/dr - (rp-1/2) (φ_i - φ_{i-1})/dr ) / dr
                # We'll implement that below directly for clarity.
                pass

            # For robustness and clarity, we'll re-implement r-term using flux formulation:
    # --- we'll rebuild using flux discretization for all interior nodes ---

def solve_poisson_pillbox(R, L, Nr, Nz, rho_func):
    """
    Solve (1/r)∂_r(r ∂φ/∂r) + ∂^2φ/∂z^2 = -rho/eps0 on 0<=r<=R, 0<=z<=L
    Dirichlet: phi=0 at r=R and z=0,L. Symmetry at r=0 (Neumann).
    Nr, Nz include boundary points.
    rho_func(r_grid, z_grid) should accept meshgrid arrays and return rho.
    Returns r, z, phi (2D arrays of shape (Nz, Nr) with phi in full grid points)
    """
    # build grids
    r = np.linspace(0, R, Nr)
    z = np.linspace(0, L, Nz)
    dr = r[1] - r[0]
    dz = z[1] - z[0]
    nr = Nr - 2
    nz = Nz - 2
    N = nr * nz

    def ij_to_k(i_r, j_z):
        return j_z * nr + i_r

    rows = []
    cols = []
    data = []
    b = np.zeros(N)

    # precompute rp+1/2 and rp-1/2 (flux form)
    for j in range(nz):
        zj = z[j+1]
        for i in range(nr):
            ir = i
            rp = r[i+1]
            k = ij_to_k(i, j)

            # flux coefficients (r-direction) using r_{i+1/2} = rp + dr/2 etc.
            rph = rp + 0.5*dr
            rpm = rp - 0.5*dr

            # coefficients multiply neighbor phi values
            # radial contribution: (1/rp) * ( r_{i+1/2}*(phi_{i+1}-phi_i)/dr - r_{i-1/2}*(phi_i - phi_{i-1})/dr ) / dr
            # => coefficients:
            a_r_plus  = rph / (rp * dr * dr)
            a_r_minus = rpm / (rp * dr * dr)
            a_r_center = -(a_r_plus + a_r_minus)

            # axial contributions (standard second difference)
            a_z_plus = 1.0 / (dz*dz)
            a_z_minus = 1.0 / (dz*dz)
            a_z_center = -(a_z_plus + a_z_minus)

            center = a_r_center + a_z_center
            rows.append(k); cols.append(k); data.append(center)

            # r+1 neighbor (i+1) -> if i == nr-1 then this neighbor is the boundary r=R (Dirichlet phi=0)
            if i < nr-1:
                rows.append(k); cols.append(ij_to_k(i+1, j)); data.append(a_r_plus)
            else:
                # neighbor is Dirichlet phi=0 -> move to RHS (nothing to add since phi_boundary=0)
                pass

            # r-1 neighbor (i-1) -> if i==0, this corresponds to r=0 symmetry point (Neumann), handled by ghost symmetry
            if i > 0:
                rows.append(k); cols.append(ij_to_k(i-1, j)); data.append(a_r_minus)
            else:
                # i == 0: neighbor is symmetry at r=0. Implement symmetry by replacing phi_{i-1} with phi_{i+1}.
                # That effectively adds a_r_minus * phi_{i+1} to RHS matrix coefficients.
                # So we move the a_r_minus coefficient to the (i+1) entry:
                if nr >= 2:
                    rows.append(k); cols.append(ij_to_k(i+1, j)); data.append(a_r_minus)
                else:
                    # degenerate small grid: skip
                    pass

            # z+1 neighbor (j+1) -> if j == nz-1 that's z=L boundary (Dirichlet -> phi=0)
            if j < nz-1:
                rows.append(k); cols.append(ij_to_k(i, j+1)); data.append(a_z_plus)
            else:
                pass

            # z-1 neighbor (j-1) -> if j==0 that's bottom boundary z=0 Dirichlet
            if j > 0:
                rows.append(k); cols.append(ij_to_k(i, j-1)); data.append(a_z_minus)
            else:
                pass

            # RHS from charge: b = -rho/eps0 evaluated at node (rp, zj)
            rho_val = rho_func(rp, zj)
            b[k] = -rho_val / eps0

    A = sp.csr_matrix((data, (rows, cols)), shape=(N, N))

    # solve
    phi_in = spla.spsolve(A, b)

    # build full phi including boundaries
    PHI = np.zeros((Nz, Nr))
    # interior
    for j in range(nz):
        for i in range(nr):
            PHI[j+1, i+1] = phi_in[ij_to_k(i, j)]
    # boundaries: phi=0 on r=R and z=0,L already zero
    # symmetry at r=0 -> set phi[:,0] = phi[:,1] (Neumann)
    PHI[:, 0] = PHI[:, 1]

    # return grids in (r,z) mesh order convenient for plotting (Z rows, R columns)
    Rgrid, Zgrid = np.meshgrid(r, z, indexing='xy')  # shapes (Nz, Nr)
    return r, z, PHI, Rgrid, Zgrid

# -----------------------
# Example usage: Gaussian charge along axis
# -----------------------
if __name__ == "__main__":
    R = 0.1   # 10 cm radius pillbox
    L = 0.05  # 5 cm length
    Nr = 101
    Nz = 81

    # Gaussian charge density centered at axis and mid-plane:
    sigma = 0.005  # radial width (m)
    sz = 0.01      # axial width (m)
    q0 = 1e-9      # total charge scaling
    def rho_func(rp, zp):
        r0 = 0.0
        z0 = L/2.0
        val = q0 * np.exp(-((rp - r0)**2)/(2*sigma**2) - ((zp - z0)**2)/(2*sz**2))
        # Normalize so total integral ~ q0 (approx)
        return val / (2*np.pi * sigma**2 * sz * np.sqrt(2*np.pi))  # not exact normalization but fine for demo

    r, z, PHI, Rg, Zg = solve_poisson_pillbox(R, L, Nr, Nz, rho_func)

    # plot potential
    plt.figure(figsize=(6,5))
    plt.pcolormesh(r, z, PHI, shading='auto')
    plt.colorbar(label='Phi (V)')
    plt.xlabel('r (m)')
    plt.ylabel('z (m)')
    plt.title('Electrostatic potential in pillbox')
    plt.show()

    # radial cut through axis at mid-plane
    mid_j = Nz//2
    plt.figure()
    plt.plot(r, PHI[mid_j, :])
    plt.xlabel('r (m)')
    plt.ylabel('Phi (V)')
    plt.title('Phi(r) at z = L/2')
    plt.grid(True)
    plt.show()

# Fixed plotting: use 1D r array for radial cut. Re-run only the plotting part using previously computed phi.
# To be safe, recompute eigenmode and produce snapshots with corrected radial cut.
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla, matplotlib.pyplot as plt
c = 299792458.0

def assemble_laplacian_cylindrical(R, L, Nr, Nz):
    r = np.linspace(0, R, Nr)
    z = np.linspace(0, L, Nz)
    dr = r[1] - r[0]; dz = z[1] - z[0]
    nr = Nr - 2; nz = Nz - 2
    N = nr * nz
    def ij_to_k(i_r, j_z): return j_z * nr + i_r
    rows, cols, data = [], [], []
    for j in range(nz):
        for i in range(nr):
            rp = r[i+1]; k = ij_to_k(i, j)
            rph = rp + 0.5*dr; rpm = rp - 0.5*dr
            a_r_plus  = rph / (rp * dr * dr); a_r_minus = rpm / (rp * dr * dr)
            a_r_center = -(a_r_plus + a_r_minus)
            a_z_plus = 1.0 / (dz*dz); a_z_minus = 1.0 / (dz*dz)
            a_z_center = -(a_z_plus + a_z_minus)
            center = a_r_center + a_z_center
            rows.append(k); cols.append(k); data.append(center)
            if i < nr-1: rows.append(k); cols.append(ij_to_k(i+1, j)); data.append(a_r_plus)
            if i > 0: rows.append(k); cols.append(ij_to_k(i-1, j)); data.append(a_r_minus)
            else:
                if nr >= 2: rows.append(k); cols.append(ij_to_k(i+1, j)); data.append(a_r_minus)
            if j < nz-1: rows.append(k); cols.append(ij_to_k(i, j+1)); data.append(a_z_plus)
            if j > 0: rows.append(k); cols.append(ij_to_k(i, j-1)); data.append(a_z_minus)
    Lmat = sp.csr_matrix((data, (rows, cols)), shape=(N, N))
    return Lmat, r, z, nr, nz, ij_to_k

def compute_eigenmode(R, L, Nr=101, Nz=81, nev=6, which_mode=0):
    Lmat, r, z, nr, nz, ij_to_k = assemble_laplacian_cylindrical(R, L, Nr, Nz)
    M = -Lmat
    eigvals, eigvecs = spla.eigsh(M, k=nev, sigma=0.0, which='LM', tol=1e-8)
    idx = np.argsort(eigvals); eigvals = eigvals[idx]; eigvecs = eigvecs[:, idx]
    phi = np.zeros((Nz, Nr)); vec = eigvecs[:, which_mode]
    for j in range(nz):
        for i in range(nr):
            phi[j+1, i+1] = vec[ij_to_k(i, j)]
    phi[:,0] = phi[:,1]
    freq = c * np.sqrt(eigvals[which_mode]) / (2*np.pi)
    Rg, Zg = np.meshgrid(r, z, indexing='xy')
    return phi, eigvals, freq, Rg, Zg, r, z

def compute_E_from_phi(phi, r, z):
    dr = r[1] - r[0]; dz = z[1] - z[0]
    dphi_dr = np.gradient(phi, dr, axis=1); dphi_dz = np.gradient(phi, dz, axis=0)
    Er = -dphi_dr; Ez = -dphi_dz; Er[:,0]=0.0
    return Er, Ez

# parameters
R = 0.05; L = 0.05; Nr = 101; Nz = 81; which_mode = 0
phi, eigvals, freq, Rg, Zg, r_arr, z_arr = compute_eigenmode(R, L, Nr, Nz, nev=6, which_mode=which_mode)
print(f"Mode {which_mode} freq = {freq:.6e} Hz")
Er_field, Ez_field = compute_E_from_phi(phi, r_arr, z_arr)
mag = np.sqrt(Er_field**2 + Ez_field**2); mag_max = mag.max()
Er_plot = Er_field / (mag_max + 1e-18); Ez_plot = Ez_field / (mag_max + 1e-18)

# Plot snapshots
phases = [0.0, np.pi/4, np.pi/2, np.pi]
fig, axes = plt.subplots(2,2, figsize=(10,8))
axes = axes.ravel()
for ax, ph in zip(axes, phases):
    scale = np.cos(ph)
    U = Er_plot * scale; V = Ez_plot * scale
    pcm = ax.pcolormesh(Rg, Zg, phi, shading='auto')
    ax.streamplot(Rg, Zg, U, V, density=1.2, linewidth=1, arrowsize=1)
    ax.set_title(f'phase = {ph:.2f} rad, scale={scale:.2f}'); ax.set_xlabel('r (m)'); ax.set_ylabel('z (m)')
    ax.set_aspect('equal', adjustable='box')
fig.colorbar(pcm, ax=axes.tolist(), orientation='vertical', fraction=0.02)
plt.tight_layout()
plt.show()

# radial cut at midplane using 1D r array
mid_j = Nz//2
plt.figure()
plt.plot(r_arr, phi[mid_j, :])
plt.xlabel('r (m)'); plt.ylabel('Ez (arb)'); plt.title('Ez(r) at z=L/2'); plt.grid(True)
plt.show()

print("Done. You can change `which_mode`, grid resolution, or cavity R/L and re-run to explore further.")
