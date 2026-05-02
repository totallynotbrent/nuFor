# Unstructured mesh research

## Question

The 2D and 3D solvers live on structured Cartesian grids: a box, a uniform
sweep, trivial indexing. Real aerospace bodies are nothing like that -- an
airfoil, a cowl, a turbine passage -- so eventually nuFor will need to solve on
unstructured meshes. The question is which representation to adopt and what it
costs to build it on the existing finite-volume core.

## Scope

This is a research note, not an implementation: it fixes the vocabulary and the
shape an unstructured mesher would take, and it sets which parts of the current
solver survive unchanged. The milestone that follows prototypes the FV
infrastructure against a small mesh.

## Findings

Three representations dominate practical CFD:

1. **Cell-centered finite volume on arbitrary polygons/polyhedra.** The mesh is
   a set of cells and the solver stores one conserved state per cell; the
   convective flux across any face uses the two adjacent cell states with the
   normal chosen by the face. Everything the solver already does -- HLLC with a
   face-normal reconstruction, MUSCL via neighbours, a CFL time step -- carries
   over essentially verbatim; what changes is the *topology* (who is adjacent to
   whom) instead of the *geometry* (a fixed grid).

2. **Vertex/median-dual schemes.** Control volumes are built around mesh nodes.
   More accurate gradients, but the dual construction is finicky to get right.

3. **Structured/blocked with a mapping.** Curvilinear grids (a structured mesh
   mapped onto a body). Keeps the index bookkeeping but adds metric terms and
   a hard meshing problem for complex bodies.

For this codebase the natural choice is **option 1, cell-centered FV on
arbitrary polyhedra**: it is the smallest change from the finite-volume core,
keeps the goal of "the same solver, a different mesh," and is the most common
representation in open-source CFD.

## What an unstructured FV implementation needs

- a fallback ordering of faces and cells and their adjacency (half-faces: an
  ordered list of the face-owning cell and the face-neighbour, or a boundary
  tag),
- face area vectors (the signed oriented area, which replaces `dx,dy,dz` as
  the metric),
- a volume per cell (the metric on the other side of the flux integral),
- a normals-into-|v|-dot-n reconstruction so the axis-aligned HLLC becomes the
  same flux with a rotated normal,
- neighbour-based gradient reconstruction for the MUSCL limiter.

The flux function itself and the conservative update are unchanged -- the axs-
aligned face is the limit of an arbitrary face with a normal pointing along an
axis, so the existing HLLC generalises by rotating the flux by the face normal.

## Decision

Adopt cell-centered finite volume on arbitrary polyhedra for the unstructured
path, reusing the HLLC flux with a face-normal rotation and a neighbour-based
limiter. Prototype it in the next milestone against a small manufactured mesh
(tetrahedra or a hand-built hexa) and verify it reproduces the structured
solver's Sod answer on a Cartesian-equivalent mesh.

## Verification implications

The unstructured path must reproduce the structured result on a mesh that is
topologically unstructured but geometrically Cartesian-aligned -- the same
sanity check the 3D solver used (a 1D problem in a box). Mass must be conserved
to round-off on a closed box, exactly as on the structured grid.

## Open questions

- Gradients for the limiter: least-squares via face neighbours is standard but
  costs a small solve per cell; is a simpler vertex-averaged gradient enough at
  first-order?
- Boundary representation: how much of the boundary-matching surface-normal
  logic lives in the mesh vs the solver?

## Date

2026-05-02