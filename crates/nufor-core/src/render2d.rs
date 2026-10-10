//! colormap rendering of a 2d scalar field to a png.

use crate::solver::ConservedState;
use crate::state::cons_to_prim;
use png::ColorType;

/// map a normalized value in [0, 1] to an rgb colour on a blue-to-red scale.
pub fn colormap(t: f64) -> [u8; 3] {
    let tt = t.clamp(0.0, 1.0);
    let stops: [(f64, [f64; 3]); 6] = [
        (0.00, [0.19, 0.08, 0.55]),
        (0.20, [0.19, 0.57, 0.85]),
        (0.40, [0.13, 0.84, 0.71]),
        (0.60, [0.57, 0.92, 0.25]),
        (0.80, [0.99, 0.66, 0.14]),
        (1.00, [0.84, 0.11, 0.11]),
    ];
    let mut rgb = [0.0f64; 3];
    for i in 0..5 {
        if tt >= stops[i].0 && tt <= stops[i + 1].0 {
            let f = (tt - stops[i].0) / (stops[i + 1].0 - stops[i].0 + 1e-12);
            let (lo, hi) = (stops[i].1, stops[i + 1].1);
            rgb = core::array::from_fn(|c| lo[c] + f * (hi[c] - lo[c]));
            break;
        }
    }
    [
        (rgb[0] * 255.0).round() as u8,
        (rgb[1] * 255.0).round() as u8,
        (rgb[2] * 255.0).round() as u8,
    ]
}

/// render an n-by-n row-major field (j increasing upward) as a png byte vector.
pub fn render_png(field: &[f64], n: usize, lo: f64, hi: f64) -> Result<Vec<u8>, String> {
    render_png_rect(field, n, n, lo, hi)
}

/// render an nx-by-ny row-major field (j increasing upward) as a png.
pub fn render_png_rect(
    field: &[f64],
    nx: usize,
    ny: usize,
    lo: f64,
    hi: f64,
) -> Result<Vec<u8>, String> {
    let range = (hi - lo).max(1e-30);
    let mut px = vec![0u8; nx * ny * 3];
    for j in 0..ny {
        // png rows run top to bottom; the field stores j = 0 at the bottom.
        let dst = ny - 1 - j;
        for i in 0..nx {
            let t = ((field[j * nx + i] - lo) / range).clamp(0.0, 1.0);
            let [r, g, b] = colormap(t);
            let base = (dst * nx + i) * 3;
            px[base] = r;
            px[base + 1] = g;
            px[base + 2] = b;
        }
    }
    let mut bytes = Vec::new();
    {
        let mut enc = png::Encoder::new(&mut bytes, nx as u32, ny as u32);
        enc.set_color(ColorType::Rgb);
        enc.set_depth(png::BitDepth::Eight);
        let mut writer = enc.write_header().map_err(|e| format!("png header: {e}"))?;
        writer
            .write_image_data(&px)
            .map_err(|e| format!("png data: {e}"))?;
    }
    Ok(bytes)
}

/// a 1d field as a line plot png (x horizontal, value vertical).
pub fn render_line1d(centers: &[f64], st: &ConservedState, gamma: f64, field: &str) -> Vec<u8> {
    let n = st.rho.len();
    let (u, _et) = cons_to_prim(&st.rho, &st.m, &st.e).unwrap_or_else(|_| {
        (
            vec![0.0; n],
            vec![st.e.iter().cloned().fold(0.0, f64::max); n],
        )
    });
    let vals: Vec<f64> = match field {
        "u" => u,
        "p" => (0..n)
            .map(|i| (gamma - 1.0) * (st.e[i] - 0.5 * st.m[i] * st.m[i] / st.rho[i].max(1e-12)))
            .collect(),
        _ => st.rho.clone(),
    };
    let w = 640usize;
    let h = 320usize;
    let bg = [16, 16, 16u8];
    let ax = [43, 43, 43u8];
    let amber = [224, 164, 88u8];
    let mut px = vec![bg[0]; w * h * 3];
    for x in 0..w {
        for c in 0..3 {
            px[(h - 1) * w * 3 + x * 3 + c] = ax[c];
            px[((h / 4) * w + x) * 3 + c] = [30, 30, 30u8][c];
            px[((h / 2) * w + x) * 3 + c] = [24, 24, 24u8][c];
            px[((3 * h / 4) * w + x) * 3 + c] = [30, 30, 30u8][c];
        }
    }
    let v0 = vals.iter().cloned().fold(f64::INFINITY, f64::min);
    let v1 = vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let span = (v1 - v0).max(1e-9);
    let x0 = *centers.first().unwrap_or(&0.0);
    let x1 = *centers.last().unwrap_or(&1.0);
    let xs = (x1 - x0).max(1e-9);
    for i in 0..n {
        let cx = (((centers[i] - x0) / xs) * (w - 2) as f64 + 1.0) as usize;
        let cy = (h as f64 - 4.0 - ((vals[i] - v0) / span) * (h - 10) as f64) as usize;
        let cy = cy.clamp(1, h - 2);
        for dx in 0..2 {
            let px_x = (cx + dx).min(w - 1);
            for c in 0..3 {
                px[cy * w * 3 + px_x * 3 + c] = amber[c];
            }
        }
    }
    let mut bytes = Vec::new();
    {
        let mut enc = png::Encoder::new(&mut bytes, w as u32, h as u32);
        enc.set_color(png::ColorType::Rgb);
        enc.set_depth(png::BitDepth::Eight);
        if let Ok(mut writer) = enc.write_header() {
            let _ = writer.write_image_data(&px);
        }
    }
    bytes
}
