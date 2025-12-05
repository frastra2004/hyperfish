import numpy as np
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