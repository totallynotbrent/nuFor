# generate the refined stagnation mesh: x clustered to 100 um at the wall x=2
# (growth 1.12 both sides), y clustered to 2 mm at the stagnation line (growth
# 1.13, mirrored to +-6 m). resolves the ~11 mm tangential layer the uniform-y
# mesh aliased. scratch tool, deleted after the case lands.

first_x, gx = 100e-6, 1.12
x0, x1, xw = 0.0, 3.0, 2.0

# upstream ramp: start AT the wall, march toward x0, then reverse so the
# list reads x0 -> xw with the small cells adjacent to the wall.
hs_up = []
x = xw
h = first_x
while x - h > x0:
    hs_up.append(h); x -= h; h *= gx
hs_up.append(x - x0)  # last (coarse) cell closes onto x0 exactly
hs_up.reverse()
faces_x = [x0]
for h in hs_up:
    faces_x.append(faces_x[-1] + h)
faces_x[-1] = xw

# downstream ramp: from the wall to x1, small at the wall, written as
# cumulative face coordinates directly.
hs_dn = []
x = xw
h = first_x
while x + h < x1:
    hs_dn.append(x + h)
    x += h
    h *= gx
hs_dn.append(x1)
faces_x = faces_x + hs_dn

# y ramp from the stagnation line (0) to 6 m: small cells at y=0.
first_y, gy = 2e-3, 1.13
hs_y = []
y = 0.0
h = first_y
while y + h < 6.0:
    hs_y.append(h); y += h; h *= gy
hs_y.append(6.0 - y)  # close onto 6.0 exactly
faces_half = [0.0]
for h in hs_y:
    faces_half.append(faces_half[-1] + h)
faces_half[-1] = 6.0
faces_y = [-v for v in reversed(faces_half)] + faces_half[1:]

with open("cases/cylinder-stagnation-fay/mesh_faces_refined.txt", "w") as f:
    for v in faces_x:
        f.write(f"{v:.10e}\n")
    f.write("\n")
    for v in faces_y:
        f.write(f"{v:.10e}\n")
print(f"x faces: {len(faces_x)} (nx={len(faces_x)-1})")
print(f"  first cell at wall (upstream): {1e6*(faces_x[1]-faces_x[0]):.0f} um" if False else
      f"  wall-adjacent cell: {1e6*(xw - faces_x[-2]):.0f} um upstream, {1e6*(faces_x[2] - xw):.0f} um downstream")
print(f"y faces: {len(faces_y)} (ny={len(faces_y)-1}), first dy at y=0: {1e3*(faces_y[1]-faces_y[0]):.1f} mm")
cells_11mm = sum(1 for a, b in zip(faces_y[:-1], faces_y[1:]) if (b-a) < 5.5e-3 and a > -0.06 and b < 0.06)
print(f"  y cells under 5.5mm near the line: {cells_11mm}")
