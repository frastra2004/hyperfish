import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

eps0 = 8.854187817e-12

def solve_poisson_pillbox(R, L, Nr, Nz, rho_func,
                          bcd_Dirichlet='auto',    # 'auto' or (phi_rR, phi_z0, phi_zL) arrays
                          bcd_Neumann='auto'):     # 'auto' or array of length Nz giving dphi/dr at r=0
    """
    Solve (1/r)∂_r(r ∂φ/∂r) + ∂^2φ/∂z^2 = -rho/eps0 on 0<=r<=R, 0<=z<=L
    Boundary condition inputs:
      - bcd_Dirichlet: 'auto' (default) => zeros at r=R and z=0,L,
                        or tuple/list (phi_rR, phi_z0, phi_zL) arrays
                      phi_rR: length Nz         (values along r=R for z[0..Nz-1])
                      phi_z0: length Nr         (values along z=0 for r[0..Nr-1])
                      phi_zL: length Nr         (values along z=L for r[0..Nr-1])
      - bcd_Neumann: 'auto' => symmetry (dphi/dr = 0) at r=0,
                     or array of length Nz giving dphi/dr at r=0 for each z.
    Nr, Nz include boundary points.
    rho_func(r, z) should accept scalars (rp, zj) or be vectorized over arrays.
    Returns r, z, PHI, Rgrid, Zgrid
    """
    # build grids
    r = np.linspace(0, R, Nr)
    z = np.linspace(0, L, Nz)
    dr = r[1] - r[0]
    dz = z[1] - z[0]
    nr = Nr - 2
    nz = Nz - 2
    N = nr * nz
    if N <= 0:
        raise ValueError("Need at least 3 points in each direction (including boundaries).")

    # Process Dirichlet BC arrays
    if bcd_Dirichlet == 'auto':
        # default zeros on right and top/bottom
        phi_rR = np.zeros(Nz, dtype=float)   # along boundary r=R sampled at all z grid points
        phi_z0 = np.zeros(Nr, dtype=float)   # along z=0 sampled at all r grid points
        phi_zL = np.zeros(Nr, dtype=float)   # along z=L sampled at all r grid points
    else:
        # expect tuple/list of three arrays
        try:
            phi_rR, phi_z0, phi_zL = bcd_Dirichlet
            phi_rR = np.asarray(phi_rR, dtype=float).ravel()
            phi_z0 = np.asarray(phi_z0, dtype=float).ravel()
            phi_zL = np.asarray(phi_zL, dtype=float).ravel()
        except Exception as e:
            raise ValueError("bcd_Dirichlet must be 'auto' or a tuple/list (phi_rR, phi_z0, phi_zL).") from e
        if phi_rR.size != Nz:
            raise ValueError(f"phi_rR must have length Nz={Nz} (got {phi_rR.size})")
        if phi_z0.size != Nr:
            raise ValueError(f"phi_z0 must have length Nr={Nr} (got {phi_z0.size})")
        if phi_zL.size != Nr:
            raise ValueError(f"phi_zL must have length Nr={Nr} (got {phi_zL.size})")

    # Process Neumann BC at r=0
    if bcd_Neumann == 'auto' or bcd_Neumann is None:
        g_r0 = np.zeros(Nz, dtype=float)   # dphi/dr at r=0 -> default 0 (symmetry)
    else:
        g_r0 = np.asarray(bcd_Neumann, dtype=float).ravel()
        if g_r0.size != Nz:
            raise ValueError(f"bcd_Neumann array must have length Nz={Nz}")

    # helper index mapping
    def ij_to_k(i_r, j_z):
        return j_z * nr + i_r

    rows = []
    cols = []
    data = []
    b = np.zeros(N, dtype=float)

    # assemble flux-form stencil and move BC contributions to RHS where needed
    for j in range(nz):
        zj = z[j+1]         # physical z for interior row j
        for i in range(nr):
            rp = r[i+1]     # physical r for interior col i
            k = ij_to_k(i, j)

            # flux coefficients in r using r_{i+1/2}, r_{i-1/2}
            rph = rp + 0.5*dr
            rpm = rp - 0.5*dr

            a_r_plus  = rph / (rp * dr * dr)
            a_r_minus = rpm / (rp * dr * dr)
            a_r_center = -(a_r_plus + a_r_minus)

            # axial second-difference
            a_z_plus = 1.0 / (dz*dz)
            a_z_minus = 1.0 / (dz*dz)
            a_z_center = -(a_z_plus + a_z_minus)

            center = a_r_center + a_z_center
            rows.append(k); cols.append(k); data.append(center)

            # r+1 neighbor
            if i < nr-1:
                rows.append(k); cols.append(ij_to_k(i+1, j)); data.append(a_r_plus)
            else:
                # neighbor is boundary at r = R: Dirichlet value phi_rR at z index j+1
                phi_b = phi_rR[j+1]   # phi at (r=R, z = zj)
                # move term a_r_plus * phi_b to RHS: subtract
                b[k] -= a_r_plus * phi_b

            # r-1 neighbor (i>0) OR r=0 treatment
            if i > 0:
                rows.append(k); cols.append(ij_to_k(i-1, j)); data.append(a_r_minus)
            else:
                # i == 0 -> uses r=0 condition.
                # Use ghost relation from Neumann at r=0:
                # (phi_1 - phi_{-1}) / (2 dr) = g_r0_at_z  => phi_{-1} = phi_1 - 2 dr * g0
                g0 = g_r0[j+1]   # Neumann at r=0 for this z
                # a_r_minus * phi_{i-1} -> a_r_minus*(phi_1 - 2 dr * g0)
                # add coefficient a_r_minus to phi_{i+1} (i+1 == 1)
                if nr >= 2:
                    rows.append(k); cols.append(ij_to_k(i+1, j)); data.append(a_r_minus)
                # move constant piece -2 dr * a_r_minus * g0 to RHS:
                b[k] += 2.0 * dr * a_r_minus * g0

            # z+1 neighbor
            if j < nz-1:
                rows.append(k); cols.append(ij_to_k(i, j+1)); data.append(a_z_plus)
            else:
                # neighbor is top boundary at z=L (Dirichlet phi_zL at r index i+1)
                phi_top = phi_zL[i+1]
                b[k] -= a_z_plus * phi_top

            # z-1 neighbor
            if j > 0:
                rows.append(k); cols.append(ij_to_k(i, j-1)); data.append(a_z_minus)
            else:
                # neighbor is bottom boundary at z=0 (Dirichlet phi_z0 at r index i+1)
                phi_bot = phi_z0[i+1]
                b[k] -= a_z_minus * phi_bot

            # RHS source term from rho
            rho_val = rho_func(rp, zj)
            b[k] += -rho_val / eps0

    A = sp.csr_matrix((data, (rows, cols)), shape=(N, N))

    # solve linear system
    phi_in = spla.spsolve(A, b)

    # build full phi including boundaries
    PHI = np.zeros((Nz, Nr), dtype=float)
    # interior
    for j in range(nz):
        for i in range(nr):
            PHI[j+1, i+1] = phi_in[ij_to_k(i, j)]

    # fill Dirichlet boundaries from provided arrays:
    # right boundary r=R -> phi_rR array length Nz (z index 0..Nz-1)
    PHI[:, -1] = phi_rR
    # bottom boundary z=0 -> phi_z0 length Nr (r index 0..Nr-1)
    PHI[0, :] = phi_z0
    # top boundary z=L -> phi_zL
    PHI[-1, :] = phi_zL

    # symmetry / Neumann at r=0:
    # If user provided Neumann g_r0 != 0, reflect ghost relation to compute phi[:,0] from phi[:,1]:
    # Use forward difference: (phi_1 - phi_(-1)) / (2 dr) = g0 and phi_(-1) = phi_1 - 2 dr g0,
    # but we need phi_0 (at r=0). A common conservative choice is to set phi_0 = phi_1 - dr * g0
    # (linear extrapolation with derivative). That gives second-order accuracy for small dr.
    # We'll use phi_0 = phi_1 - dr * g0.
    for jj in range(Nz):
        g0_full = g_r0[jj]   # Neumann at r=0 for this z
        PHI[jj, 0] = PHI[jj, 1] - dr * g0_full

    # generate grid meshes for output
    Rgrid, Zgrid = np.meshgrid(r, z, indexing='xy')  # shapes (Nz, Nr)
    return r, z, PHI, Rgrid, Zgrid
