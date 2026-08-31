"""Shared FEA helpers for the Magneto IEM mechanical validation.

Everything geometric is pulled straight out of ``generate.py`` -- the same SDF
functions that write the STLs -- so the FEA and the printed part cannot drift
apart.  Nothing here re-implements the geometry.

Meshing strategy
----------------
Two meshers are provided:

``voxel_hex``    Structured trilinear-hex ("voxel") mesh: sample the SDF on a
                 regular grid, keep the elements whose centroid is inside.
                 This is the standard approach for TPMS / trabecular-bone
                 micro-FE.  It is robust (no surface-meshing failures on a
                 0.2 mm gyroid wall) and it is known to be mildly *stiff*-biased
                 and to produce spurious stress spikes at the staircase corners
                 -- both are reported honestly at the call sites.

``tet_from_sdf`` Marching cubes -> STL -> gmsh volume mesh.  Used for the smooth
                 macro envelopes (wing, jacket skin) where the surface meshes
                 cleanly and stress recovery is better behaved.

Material: Ti-6Al-4V, E = 110 GPa, nu = 0.31, sigma_y = 900 MPa (LPBF as-printed
is usually quoted 950-1100 MPa; 900 is the conservative wrought-anneal figure).
Units throughout: mm, N, MPa (= N/mm^2).
"""

import copy
import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))

import generate as gen  # noqa: E402  -- the real generator, not a copy

from skfem import (Basis, ElementHex1, ElementTetP1, ElementVector,  # noqa: E402
                   MeshHex, MeshTet, asm)
from skfem.helpers import sym_grad  # noqa: E402
from skfem.models.elasticity import lame_parameters, linear_elasticity  # noqa: E402

# ---------------------------------------------------------------------------
# material
# ---------------------------------------------------------------------------

E_TI = 110.0e3      # MPa
NU_TI = 0.31
SY_TI = 900.0       # MPa, yield
SY_HALF = 450.0     # MPa, "fatigue-relevant" flag level asked for
RHO_TI = 4.43e-3    # g/mm^3


def geom(**overrides):
    """A fresh generate.G built from the shipped PARAMS."""
    P = copy.deepcopy(gen.PARAMS)
    P.update(overrides)
    return gen.G(P), P


# ---------------------------------------------------------------------------
# meshing
# ---------------------------------------------------------------------------

def sample_sdf(fn, bounds, h):
    """Sample an SDF on a regular grid.  Returns (field, x, y, z)."""
    (x0, x1), (y0, y1), (z0, z1) = bounds
    x = np.arange(x0, x1 + 0.5 * h, h)
    y = np.arange(y0, y1 + 0.5 * h, h)
    z = np.arange(z0, z1 + 0.5 * h, h)
    C = gen.Ctx(x, y, z)
    f = fn(x.reshape(-1, 1, 1), y.reshape(1, -1, 1), z.reshape(1, 1, -1), C)
    return np.asarray(f, dtype=np.float64), x, y, z


def voxel_hex(fn, bounds, h, level=0.0):
    """Structured hex mesh of {fn < level}, elements picked by centroid.

    Returns (MeshHex, centroids (3, nelem)) with the largest connected component
    only -- stray islands from a marching lattice would make K singular.
    """
    (x0, x1), (y0, y1), (z0, z1) = bounds
    nx = int(round((x1 - x0) / h))
    ny = int(round((y1 - y0) / h))
    nz = int(round((z1 - z0) / h))
    # element centres
    cx = x0 + h * (np.arange(nx) + 0.5)
    cy = y0 + h * (np.arange(ny) + 0.5)
    cz = z0 + h * (np.arange(nz) + 0.5)
    C = gen.Ctx(cx, cy, cz)
    f = fn(cx.reshape(-1, 1, 1), cy.reshape(1, -1, 1), cz.reshape(1, 1, -1), C)
    occ = np.asarray(f, dtype=np.float64) < level
    occ = _largest_component(occ)

    m = MeshHex.init_tensor(x0 + h * np.arange(nx + 1),
                            y0 + h * np.arange(ny + 1),
                            z0 + h * np.arange(nz + 1))
    # init_tensor's element order is NOT C-order, so map back through centroids
    c = m.p[:, m.t].mean(axis=1)
    ii = np.clip(np.round((c[0] - x0) / h - 0.5).astype(int), 0, nx - 1)
    jj = np.clip(np.round((c[1] - y0) / h - 0.5).astype(int), 0, ny - 1)
    kk = np.clip(np.round((c[2] - z0) / h - 0.5).astype(int), 0, nz - 1)
    keep = np.where(occ[ii, jj, kk])[0]
    if keep.size == 0:
        raise RuntimeError("voxel_hex: nothing inside the level set")
    m = m.restrict(keep)
    cent = m.p[:, m.t].mean(axis=1)
    return m, cent


def _largest_component(occ):
    """Largest face-connected (6-neighbour) component of a boolean voxel array.

    Face connectivity, not 26, so the retained mesh cannot hang together through
    a single shared node or edge -- that would be a mechanism and would make K
    singular.
    """
    from scipy import ndimage
    lab, n = ndimage.label(occ)
    if n <= 1:
        return occ
    sizes = ndimage.sum(occ, lab, range(1, n + 1))
    best = int(np.argmax(sizes)) + 1
    dropped = occ.sum() - sizes[best - 1]
    if dropped:
        print(f"    (dropped {int(dropped)} disconnected voxels of {int(occ.sum())})")
    return lab == best


def tet_from_sdf(fn, bounds, h_mc, h_tet, tag="part"):
    """Marching cubes on the SDF, then a gmsh tetrahedral volume mesh."""
    import gmsh
    import trimesh

    f, x, y, z = sample_sdf(fn, bounds, h_mc)
    f = np.pad(f, 1, mode="constant", constant_values=abs(h_mc) * 4 + 1.0)
    from skimage import measure
    verts, faces, _, _ = measure.marching_cubes(f, level=0.0, spacing=(h_mc,) * 3)
    verts += np.array([x[0] - h_mc, y[0] - h_mc, z[0] - h_mc])
    surf = trimesh.Trimesh(vertices=verts, faces=faces, process=False)
    if surf.volume < 0:
        surf.invert()
    stl = os.path.join(_HERE, f"_tmp_{tag}.stl")
    surf.export(stl)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.merge(stl)
    gmsh.model.mesh.classifySurfaces(math.radians(40.0), True, True, math.radians(180.0))
    gmsh.model.mesh.createGeometry()
    surfs = [t for (d, t) in gmsh.model.getEntities(2)]
    loop = gmsh.model.geo.addSurfaceLoop(surfs)
    gmsh.model.geo.addVolume([loop])
    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMin", h_tet * 0.6)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h_tet)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.model.mesh.generate(3)
    ntag, ncoord, _ = gmsh.model.mesh.getNodes()
    p = np.asarray(ncoord).reshape(-1, 3).T
    order = np.argsort(ntag)
    remap = np.zeros(int(ntag.max()) + 1, dtype=int)
    remap[np.asarray(ntag, dtype=int)[order]] = np.arange(len(ntag))
    p = p[:, order]
    et, en = gmsh.model.mesh.getElementsByType(4)[0:2]
    t = remap[np.asarray(en, dtype=int).reshape(-1, 4)].T
    gmsh.finalize()
    os.remove(stl)

    used = np.unique(t)
    rem2 = -np.ones(p.shape[1], dtype=int)
    rem2[used] = np.arange(used.size)
    return MeshTet(np.ascontiguousarray(p[:, used]),
                   np.ascontiguousarray(rem2[t])), surf.volume


# ---------------------------------------------------------------------------
# linear elasticity with an element-wise (graded) modulus
# ---------------------------------------------------------------------------

def rotate_mesh(m, R):
    """Rigid-rotate a mesh by R.  Isotropic elasticity is frame-indifferent, so
    solving on the rotated mesh is exact -- and it makes an arbitrary press
    direction into the +z coordinate dof, which turns each frictionless contact
    constraint into a single Dirichlet dof instead of a multipoint constraint."""
    cls = type(m)
    return cls(np.ascontiguousarray(R @ m.p), np.ascontiguousarray(m.t))


def make_basis(m, intorder=2):
    """intorder=2 is FULL (2x2x2 Gauss) integration for a trilinear hex -- not
    reduced, so there is no hourglassing risk; it just avoids skfem's default
    4x4x4 rule, which costs 8x the assembly time for no accuracy."""
    el = ElementHex1() if isinstance(m, MeshHex) else ElementTetP1()
    return Basis(m, ElementVector(el), intorder=intorder)


def stiffness(basis, E, nu=NU_TI):
    """Assemble K.  E may be a scalar or an (nelem,) array (graded material)."""
    if np.isscalar(E):
        lam, mu = lame_parameters(float(E), nu)
        return asm(linear_elasticity(lam, mu), basis)
    E = np.asarray(E, dtype=float)
    nq = basis.X.shape[1] if hasattr(basis, "X") else 1
    lam_e = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu_e = E / (2 * (1 + nu))
    lam = np.tile(lam_e.reshape(-1, 1), (1, nq))
    mu = np.tile(mu_e.reshape(-1, 1), (1, nq))
    return asm(linear_elasticity(lam, mu), basis)


def von_mises(basis, u, E, nu=NU_TI):
    """Element-wise von Mises stress at quadrature points -> (nelem, nqp)."""
    if np.isscalar(E):
        E = np.full(basis.mesh.t.shape[1], float(E))
    eps = sym_grad(basis.interpolate(u))
    tr = eps[0, 0] + eps[1, 1] + eps[2, 2]
    lam = (E * nu / ((1 + nu) * (1 - 2 * nu))).reshape(-1, 1)
    mu = (E / (2 * (1 + nu))).reshape(-1, 1)
    s = np.empty_like(eps)
    for i in range(3):
        for j in range(3):
            s[i, j] = 2.0 * mu * eps[i, j] + (lam * tr if i == j else 0.0)
    d = s.copy()
    p = (s[0, 0] + s[1, 1] + s[2, 2]) / 3.0
    for i in range(3):
        d[i, i] = s[i, i] - p
    j2 = 0.5 * sum(d[i, j] ** 2 for i in range(3) for j in range(3))
    return np.sqrt(3.0 * j2)


def solve_fixed(K, basis, D, x_pre, tol=1e-11, verbose=False,
                direct_max=260000):
    """Solve K u = 0 with Dirichlet dofs D held at x_pre[D].

    Sparse direct (SuperLU) up to `direct_max` dof, AMG-preconditioned CG above
    it.  Direct is preferred wherever it fits: smoothed-aggregation AMG with a
    rigid-body near-nullspace converges fine on chunky solids but STALLS on thin
    shells (the 0.2 mm wing sheet), where the near-nullspace is dominated by
    bending modes that rigid-body vectors do not span.  SuperLU on the same
    system is both robust and fast (17 s at 100k dof), because a thin shell has
    a favourable elimination tree.
    Returns (u, r) where r = K u are the nodal reactions.
    """
    import scipy.sparse as sp

    n = K.shape[0]
    D = np.unique(np.asarray(D, dtype=np.int64))
    x = np.zeros(n)
    x[D] = x_pre[D]
    free = np.setdiff1d(np.arange(n), D)
    A = K[free][:, free].tocsr()
    b = -np.asarray(K @ x)[free]

    if A.shape[0] < direct_max:
        u_f = sp.linalg.splu(A.tocsc()).solve(b)
    else:
        import pyamg
        nd = basis.nodal_dofs                       # (3, nnodes)
        p = basis.mesh.p
        B = np.zeros((n, 6))
        for c in range(3):
            B[nd[c], c] = 1.0
        B[nd[0], 3] = -p[1]; B[nd[1], 3] = p[0]     # Rz
        B[nd[0], 4] = -p[2]; B[nd[2], 4] = p[0]     # Ry
        B[nd[1], 5] = -p[2]; B[nd[2], 5] = p[1]     # Rx
        # Jacobi prolongation smoothing + the default strength measure: ~200x
        # cheaper to set up than energy-minimisation smoothing here, and with
        # the rigid-body near-nullspace supplied it still converges in tens of
        # iterations on these systems.
        ml = pyamg.smoothed_aggregation_solver(A, B=B[free], smooth="jacobi",
                                               max_coarse=1000)
        res = []
        u_f = ml.solve(b, tol=tol, accel="cg", maxiter=800, residuals=res)
        rel = res[-1] / max(res[0], 1e-300)
        if verbose:
            print(f"      AMG-CG {len(res)-1} its, rel resid {rel:.2e}")
        if rel > 1e-7:
            raise RuntimeError(f"AMG-CG did not converge (rel resid {rel:.2e})")

    x[free] = u_f
    return x, np.asarray(K @ x)
