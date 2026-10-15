# generate the nose-only stagnation mesh: x clustered toward the body wall
# at x=2.0 from both sides, y uniform across the body width. the domain ends
# at the sphere's widest point (x=3.0) so the wake never exists: the
# shoulder region is a clean supersonic outflow and the frozen-floor wake
# pocket that poisoned the full-wake march never forms.
# scratch tool, deleted after the case lands.
import math

ny = 384
x0, x1 = 0.0, 3.0
y0, y1 = -6.0, 6.0
xw = 2.0

# upstream side: natural geometric fill from the wall (x=2) back to the
# inflow (x=0), first cell 250 um, growth 1.18; last cell clipped to the
# remaining span so the inflow face lands exactly at 0.
first_cell = 250e-6
growth = 1.18
hs_up = []
remaining = xw - x0
h = first_cell
while remaining > h:
    hs_up.append(h)
    remaining -= h
    h *= growth
if remaining > 1e-9 and hs_up:
    hs_up.append(remaining)
hs_up.reverse()  # now ordered large-at-inflow -> small-at-wall
faces_x = [x0]
for h in hs_up:
    faces_x.append(faces_x[-1] + h)
faces_x[-1] = xw

# downstream side: the same ramp from the wall out to x=3.0 (the sphere's
# widest point), growth 1.18.
faces_x_dn = [xw]
h = first_cell
while True:
    nxt = faces_x_dn[-1] + h
    if nxt >= x1:
        faces_x_dn.append(x1)
        break
    faces_x_dn.append(nxt)
    h *= growth
faces_x = faces_x + faces_x_dn[1:]

faces_y = [y0 + (y1 - y0) * j / ny for j in range(ny + 1)]

with open("cases/cylinder-stagnation-fay/mesh_faces.txt", "w") as f:
    for v in faces_x:
        f.write(f"{v:.10e}\n")
    f.write("\n")
    for v in faces_y:
        f.write(f"{v:.10e}\n")
print(f"x faces: {len(faces_x)} (nx = {len(faces_x)-1}), "
      f"first/last x: {faces_x[0]:.3f}/{faces_x[-1]:.2f}, wall at {xw}")
print(f"first cell upstream: {hs_up[0]:.3e} m, downstream: {faces_x_dn[1]-faces_x_dn[0]:.3e}")
print(f"y faces: {len(faces_y)} (ny = {ny}), span {faces_y[-1]-faces_y[0]}")
