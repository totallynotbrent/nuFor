# liou's cdv nozzle (aiaa 87-0355): cosine area law, axisymmetric radius.
# the solid wall polygon: contour from inlet to exit, closed over the
# domain top. run from the case dir; writes mesh_faces.txt is not needed
# (the body polygon is in case.toml) but the contour table is written for
# the post-processing gates.
import math

def area(x):
    if x < 5.0:
        return 1.75 - 0.75 * math.cos((0.2 * x - 1.0) * math.pi)
    return 1.25 - 0.25 * math.cos((0.2 * x - 1.0) * math.pi)

def radius(x):
    return math.sqrt(area(x) / math.pi)

rows = []
for i in range(0, 201):
    x = 10.0 * i / 200
    rows.append((x, radius(x)))
with open("contour.csv", "w") as f:
    f.write("# x, r(x) for liou's cdv nozzle (aiaa 87-0355)\n")
    for x, r in rows:
        f.write(f"{x:.4f},{r:.6f}\n")

# the polygon verts (coarser than the table: 61 points, one every 1/6 x-unit
# plus the closure over the top):
verts = []
for i in range(0, 61):
    x = 10.0 * i / 60
    verts.append((x, radius(x)))
r_max = max(r for _, r in verts)
y_top = r_max * 1.35
verts.append((10.0, y_top))
verts.append((0.0, y_top))
with open("verts.txt", "w") as f:
    f.write("[[" + ", ".join(f"[{x:.4f}, {r:.6f}]" for x, r in verts) + "]]\n")
print(f"contour: A(0)={area(0.0):.3f} A(5)={area(5.0):.3f} A(10)={area(10.0):.3f}")
print(f"r_in={radius(0.0):.4f} r_throat={radius(5.0):.4f} r_exit={radius(10.0):.4f}")
print(f"verts: {len(verts)}")
