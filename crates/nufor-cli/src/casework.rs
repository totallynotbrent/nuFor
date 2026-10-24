//! the reader and renderer the case viewer shares with the results
//! ui: recover a finished run's state from its vtk, then paint
//! fields with the body silhouette in the house colormap.

use nufor_core::{render_png_gray, render_png_rect, ConservedState2d, Grid2d};

pub fn case_vtk_read(
    path: &std::path::Path,
    model: nufor_core::ThermoModel,
) -> Result<CaseFields, String> {
    let f = std::fs::read_to_string(path).map_err(|e| format!("read vtk: {e}"))?;
    let (nx, ny) = f
        .lines()
        .find_map(|l| {
            let l = l.trim();
            if let Some(rest) = l.strip_prefix("DIMENSIONS") {
                let mut it = rest.split_whitespace();
                let nx: usize = it.next()?.parse().ok()?;
                let ny: usize = it.next()?.parse().ok()?;
                return Some((nx, ny));
            }
            None
        })
        .ok_or_else(|| "read vtk: no DIMENSIONS".to_string())?;
    let field = |name: &str| -> Option<Vec<f64>> {
        let i = f.find(&format!("SCALARS {name}"))?;
        let lt = f[i..].find("LOOKUP_TABLE")? + i;
        let j = f[lt..].find('\n')? + lt + 1;
        Some(
            f[j..]
                .split_whitespace()
                .take(nx * ny)
                .filter_map(|x| x.parse::<f64>().ok())
                .collect(),
        )
    };
    let n = nx * ny;
    let rho = field("density").ok_or("read vtk: no density")?;
    if rho.len() != n {
        return Err("read vtk: short density block".into());
    }
    let e = field("energy").ok_or("read vtk: no energy")?;
    let mach = field("mach").unwrap_or_else(|| vec![0.0; n]);
    // the vtk carries the closure's own pressure and temperature, which
    // the species lanes cannot recompute from (rho, e) alone: the
    // composition is not in the conserved state.
    let pressure = field("pressure");
    let temperature = field("temperature");
    // recover the speed per cell: |v| = mach * a, with a from the
    // closure at the running estimate of the internal energy.
    let mut speed = vec![0.0f64; n];
    for k in 0..n {
        let et = e[k] / rho[k].max(1e-12);
        let mut ke = 0.0;
        for _ in 0..3 {
            let e_int = (et - ke).max(1.0e4);
            let a = match nufor_core::thermo_sound(model, &[rho[k]], &[e_int], &[0.0], &[0.0]) {
                Ok(v) => v[0],
                Err(_) => 300.0,
            };
            let s_k = mach[k] * a;
            ke = 0.5 * s_k * s_k;
        }
        // the iteration can overshoot in near-vacuum wake cells where
        // the closure's sound speed estimate is unstable; never let
        // the recovered kinetic energy exceed half the total.
        ke = ke.min(et * 0.5);
        speed[k] = (2.0 * ke).sqrt();
    }
    let mx: Vec<f64> = (0..n).map(|k| rho[k] * speed[k]).collect();
    let st = ConservedState2d {
        rho,
        mx,
        my: vec![0.0; n],
        e,
    };
    Ok(CaseFields {
        st,
        pressure,
        temperature,
    })
}

/// the fields a case vtk carries: the conserved state the renderer
/// works from, plus the closure's own pressure and temperature when
/// the vtk has them. the species lanes cannot recompute either from
/// (rho, e) alone, because the composition is not in the conserved
/// state.
pub struct CaseFields {
    pub st: ConservedState2d,
    pub pressure: Option<Vec<f64>>,
    pub temperature: Option<Vec<f64>>,
}

/// the viewer's variant of the render: same field math and body
/// silhouette, but the color window comes from the finished field
/// itself (2nd/98th percentile over live cells) so any case reads
/// well, and temperature and schlieren are first-class fields.
/// returns the png bytes and the color window for the colorbar.
pub fn render_with_body_ranged(
    g: &Grid2d,
    fields: &CaseFields,
    field: &str,
    body: &Option<crate::CaseBodyFns>,
    model: nufor_core::ThermoModel,
) -> (Vec<u8>, f64, f64) {
    let st = &fields.st;
    let n = g.nx * g.ny;
    let u: Vec<f64> = (0..n).map(|k| st.mx[k] / st.rho[k]).collect();
    let v: Vec<f64> = (0..n).map(|k| st.my[k] / st.rho[k]).collect();
    let e_int: Vec<f64> = (0..n)
        .map(|k| st.e[k] / st.rho[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]))
        .collect();
    // the closure's own pressure comes from the vtk when it has one:
    // the species lanes cannot recompute it from (rho, e) alone.
    let p = match fields.pressure.as_ref() {
        Some(v) if v.len() == n => v.clone(),
        _ => nufor_core::thermo_pressure(model, &st.rho, &e_int, &u, &v)
            .unwrap_or_else(|_| vec![1e-12; n]),
    };
    let p: Vec<f64> = p.iter().map(|x| x.max(1e-12)).collect();
    let a =
        nufor_core::thermo_sound(model, &st.rho, &e_int, &u, &v).unwrap_or_else(|_| vec![1.0; n]);
    let a: Vec<f64> = a.iter().map(|x| x.max(1e-6)).collect();
    let mut data: Vec<f64> = match field {
        "mach" => (0..n)
            .map(|k| (u[k] * u[k] + v[k] * v[k]).sqrt() / a[k])
            .collect(),
        "schlieren" => {
            // the classic density-gradient magnitude: |grad rho| per
            // cell (central differences, one-sided at the borders),
            // log-scaled so the shock structure reads. the caller
            // percentiles the window, so no fixed range here.
            let mut s = vec![0.0f64; n];
            for j in 0..g.ny {
                for i in 0..g.nx {
                    let k = j * g.nx + i;
                    let gx = if i > 0 && i < g.nx - 1 {
                        (st.rho[k + 1] - st.rho[k - 1]) / (2.0 * g.dx)
                    } else {
                        0.0
                    };
                    let gy = if j > 0 && j < g.ny - 1 {
                        (st.rho[(j + 1) * g.nx + i] - st.rho[(j - 1) * g.nx + i]) / (2.0 * g.dy)
                    } else {
                        0.0
                    };
                    s[k] = (gx * gx + gy * gy).sqrt().ln_1p();
                }
            }
            s
        }
        "p" => p.clone(),
        "t" => match fields.temperature.as_ref() {
            Some(v) if v.len() == n => v.clone(),
            _ => {
                let t = nufor_core::thermo_temperature(model, &st.rho, &p);
                match t {
                    Ok(v) if v.len() == n => v,
                    _ => st.rho.clone(),
                }
            }
        },
        _ => st.rho.clone(),
    };
    if data.len() != n {
        // the renderer indexes the field by cell, so a short field is a
        // wiring bug, not a rendering choice.
        return (Vec::new(), 0.0, 1.0);
    }
    // which cells are inside the body (masked from the range math)
    let mut solid = vec![false; n];
    if let Some((dist, _)) = body {
        for j in 0..g.ny {
            for i in 0..g.nx {
                let (x, y) = (g.centers_x[j * g.nx + i], g.centers_y[j * g.nx + i]);
                if dist(x, y) < 0.0 {
                    solid[j * g.nx + i] = true;
                }
            }
        }
    }
    // the color window: percentiles over live cells, with a fallback
    // when everything is flat (a uniform field renders as the low
    // color rather than dividing by zero).
    let mut live: Vec<f64> = data
        .iter()
        .zip(&solid)
        .filter(|(_, s)| !**s)
        .map(|(x, _)| *x)
        .filter(|x| x.is_finite())
        .collect();
    live.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let (lo, mut hi) = if live.is_empty() {
        (0.0, 1.0)
    } else {
        let q = |f: f64| live[((live.len() as f64 - 1.0) * f).round() as usize];
        (q(0.02), q(0.98))
    };
    if hi - lo < 1e-12 {
        hi = lo + 1.0;
    }
    // mark solid cells below the window so they paint dark
    for k in 0..n {
        if solid[k] {
            data[k] = lo - 1.0;
        }
    }
    let png = if field == "schlieren" {
        // the classic schlieren read is grayscale on a dark field:
        // black freestream, white shock bands, same |grad rho| window.
        render_png_gray(&data, g.nx, g.ny, lo, hi).unwrap_or_default()
    } else {
        render_png_rect(&data, g.nx, g.ny, lo, hi).unwrap_or_default()
    };
    (png, lo, hi)
}
