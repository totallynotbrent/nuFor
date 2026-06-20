# 3D HLLC solver

With the 3D state and grid in place, [32] adds the solver that actually moves
it: the three-wave HLLC flux on all three face families, MUSCL reconstruction,
and a two-stage time step. The validation confirms the three-dimensional path
reproduces the known one-dimensional physics.

## The flux

The 3D HLLC is the 2D scheme with a third velocity. Across any axis-aligned
face the flow is rotated into the face frame (a normal component and two
passive tangents), the two-wave star state is built the same way, and the extra
tangent momentum is carried through unchanged. `hllc_flux3` takes an `axis`
argument (0 = normal x, 1 = y, 2 = z) and returns a five-component flux. On
uniform flow it reduces to the exact physical flux -- the test checks that
directly.

## The step

`advance3d` walks the domain in three perpendicular sweeps -- one per face
family -- reconstructing each primitive strip with the van Leer limiter,
computing its HLLC fluxes, and adding the divergence to the conservative
update. `advance3d_rk2` wraps the two-stage Heun step so MUSCL's second-order
accuracy is not masked by first-order time integration.

## Verification

The Sod tube, a pure 1D problem, is run in a 3D box whose y and z directions
are uniform. If the three-dimensional machinery is right, the x-profile must be
exactly the classic Sod tube:

- the star plateau (between the rarefaction tail and the contact) reads 0.425
  against the known 0.426,
- the far-left and far-right states hold at rho = 1.0 and 0.125,
- density stays positive throughout.

All three hold, so the 3D step is a faithful extension of the verified 1D and
2D paths rather than a guess.

## See also

- [[numerics/3d|3D foundations]]
- [[numerics/hll-vs-hllc|HLL vs HLLC]]
- [[operations/memory-budget|Memory budget]]
- [[numerics/sod|Sod verification]]
