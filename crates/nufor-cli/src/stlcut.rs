//! stl import: slice a stereo-lithography triangle mesh at a plane
//! and emit the cross-section as a closed 2d polygon, so an exported
//! CAD body can drive the case file's polygon section.

/// one stl triangle: three vertices in mesh coordinates.
#[derive(Clone, Copy)]
pub struct StlTri {
    pub v: [(f64, f64, f64); 3],
}

/// parse an ascii or binary stl file into triangles.
pub fn load_stl(path: &str) -> Result<Vec<StlTri>, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("could not read {path}: {e}"))?;
    // binary stl: 80-byte header, uint32 count, 50 bytes per triangle.
    // ascii stl starts with "solid" and contains the word "facet".
    let head = &bytes[..bytes.len().min(200)];
    let is_ascii = head.starts_with(b"solid") && bytes.windows(5).take(4096).any(|w| w == b"facet");
    if is_ascii {
        parse_ascii(std::str::from_utf8(&bytes).map_err(|e| format!("stl not utf-8: {e}"))?)
    } else {
        parse_binary(&bytes)
    }
}

fn parse_ascii(text: &str) -> Result<Vec<StlTri>, String> {
    let mut tris = Vec::new();
    let mut verts: Vec<(f64, f64, f64)> = Vec::new();
    for line in text.lines() {
        let t = line.trim();
        if let Some(rest) = t.strip_prefix("vertex") {
            let mut it = rest.split_whitespace();
            let (Some(x), Some(y), Some(z)) = (it.next(), it.next(), it.next()) else {
                return Err("stl: malformed vertex line".into());
            };
            let (x, y, z): (f64, f64, f64) = (
                x.parse().map_err(|_| "stl: bad vertex number")?,
                y.parse().map_err(|_| "stl: bad vertex number")?,
                z.parse().map_err(|_| "stl: bad vertex number")?,
            );
            verts.push((x, y, z));
            if verts.len() == 3 {
                tris.push(StlTri {
                    v: [verts[0], verts[1], verts[2]],
                });
                verts.clear();
            }
        }
    }
    if tris.is_empty() {
        return Err("stl: no triangles found (ascii)".into());
    }
    Ok(tris)
}

fn parse_binary(b: &[u8]) -> Result<Vec<StlTri>, String> {
    if b.len() < 84 {
        return Err("stl: file too short for binary stl".into());
    }
    let n = u32::from_le_bytes([b[80], b[81], b[82], b[83]]) as usize;
    if b.len() < 84 + n * 50 {
        return Err("stl: truncated binary stl".into());
    }
    let mut tris = Vec::with_capacity(n);
    for i in 0..n {
        let off = 84 + i * 50;
        let rd = |k: usize| {
            f32::from_le_bytes([b[off + k], b[off + k + 1], b[off + k + 2], b[off + k + 3]])
        };
        let f = |k: usize| (rd(k) as f64, rd(k + 4) as f64, rd(k + 8) as f64);
        tris.push(StlTri {
            v: [f(12), f(24), f(36)],
        });
    }
    Ok(tris)
}

/// slice the mesh at plane z = z0 and chain the crossing segments into
/// one closed loop of 2d points. the loop is emitted starting from its
/// leftmost point, ordered counterclockwise in the (x, y) plane.
pub fn slice_at(tris: &[StlTri], z0: f64) -> Result<Vec<(f64, f64)>, String> {
    let tol2 = 1e-18f64;
    let mut segs: Vec<((f64, f64), (f64, f64))> = Vec::new();
    for t in tris {
        let d = [t.v[0].2 - z0, t.v[1].2 - z0, t.v[2].2 - z0];
        // crossing edges: endpoints on opposite sides (or one touching)
        let mut pts: Vec<(f64, f64)> = Vec::new();
        for (a, b) in [(0usize, 1usize), (1, 2), (2, 0)] {
            let (da, db) = (d[a], d[b]);
            if (da > 0.0 && db < 0.0) || (da < 0.0 && db > 0.0) {
                let s = da / (da - db);
                let p = (
                    t.v[a].0 + s * (t.v[b].0 - t.v[a].0),
                    t.v[a].1 + s * (t.v[b].1 - t.v[a].1),
                );
                pts.push(p);
            }
        }
        if pts.len() == 2 {
            let dx = pts[1].0 - pts[0].0;
            let dy = pts[1].1 - pts[0].1;
            if dx * dx + dy * dy > tol2 {
                segs.push((pts[0], pts[1]));
            }
        }
    }
    if segs.is_empty() {
        return Err(format!("stl-cut: no triangles cross z = {z0}"));
    }
    // chain segments into a loop: cluster endpoints within a
    // tolerance (stl coords are f32-mangled, exact equality is
    // unreliable), then walk the chain from the longest segment so
    // near-degenerate micro-segments get skipped.
    let tol = 1e-6;
    let d2 = |a: (f64, f64), b: (f64, f64)| {
        let (dx, dy) = (a.0 - b.0, a.1 - b.1);
        dx * dx + dy * dy
    };
    let mut coords: Vec<(f64, f64)> = Vec::new();
    for (s0, s1) in &segs {
        coords.push(*s0);
        coords.push(*s1);
    }
    let mut cluster = vec![usize::MAX; coords.len()];
    let mut ncluster = 0usize;
    for i in 0..coords.len() {
        if cluster[i] != usize::MAX {
            continue;
        }
        cluster[i] = ncluster;
        loop {
            let mut grew = false;
            for j in 0..coords.len() {
                if cluster[j] == usize::MAX {
                    for k in 0..coords.len() {
                        if cluster[k] == ncluster && d2(coords[j], coords[k]) < tol * tol {
                            cluster[j] = ncluster;
                            grew = true;
                            break;
                        }
                    }
                }
            }
            if !grew {
                break;
            }
        }
        ncluster += 1;
    }
    let mut cen = vec![(0.0f64, 0.0f64); ncluster];
    let mut cnt = vec![0usize; ncluster];
    for (i, c) in cluster.iter().enumerate() {
        cen[*c].0 += coords[i].0;
        cen[*c].1 += coords[i].1;
        cnt[*c] += 1;
    }
    for c in 0..ncluster {
        cen[c] = (cen[c].0 / cnt[c] as f64, cen[c].1 / cnt[c] as f64);
    }
    let key: Vec<(usize, usize, f64)> = segs
        .iter()
        .enumerate()
        .map(|(i, (s0, s1))| (cluster[2 * i], cluster[2 * i + 1], d2(*s0, *s1)))
        .collect();
    let mut order: Vec<usize> = (0..segs.len()).collect();
    order.sort_by(|a, b| key[*b].2.partial_cmp(&key[*a].2).unwrap());
    let start_seg = order[0];
    let mut used = vec![false; segs.len()];
    used[start_seg] = true;
    let mut loop_pts: Vec<(f64, f64)> = vec![cen[key[start_seg].0], cen[key[start_seg].1]];
    loop {
        let tail = loop_pts.last().copied().unwrap_or((0.0, 0.0));
        // the tail's cluster: nearest cluster center
        let mut tail_cluster = 0usize;
        let mut bd = f64::INFINITY;
        for (c, ctr) in cen.iter().enumerate() {
            let dd = d2(tail, *ctr);
            if dd < bd {
                bd = dd;
                tail_cluster = c;
            }
        }
        let mut found = false;
        for &i in &order {
            if used[i] {
                continue;
            }
            let (cp, cq, _) = key[i];
            if cp == tail_cluster {
                loop_pts.push(cen[cq]);
                used[i] = true;
                found = true;
                break;
            }
            if cq == tail_cluster {
                loop_pts.push(cen[cp]);
                used[i] = true;
                found = true;
                break;
            }
        }
        if !found {
            break;
        }
        if d2(*loop_pts.last().unwrap(), loop_pts[0]) < tol * tol {
            loop_pts.pop();
            break;
        }
    }
    if loop_pts.len() < 3 {
        return Err("stl-cut: the cross-section did not close into a loop".into());
    }

    // orient counterclockwise and start at the leftmost vertex
    let area: f64 = loop_pts
        .windows(2)
        .map(|w| w[0].0 * w[1].1 - w[1].0 * w[0].1)
        .sum::<f64>()
        + loop_pts[loop_pts.len() - 1].0 * loop_pts[0].1
        - loop_pts[0].0 * loop_pts[loop_pts.len() - 1].1;
    if area < 0.0 {
        loop_pts.reverse();
    }
    let start = loop_pts
        .iter()
        .enumerate()
        .min_by(|a, b| a.1 .0.partial_cmp(&b.1 .0).unwrap())
        .map(|(i, _)| i)
        .unwrap_or(0);
    loop_pts.rotate_left(start);
    Ok(loop_pts)
}
