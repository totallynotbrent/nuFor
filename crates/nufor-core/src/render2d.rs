//! colormap rendering of a 2d scalar field to a png.

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
    let range = (hi - lo).max(1e-30);
    let mut px = vec![0u8; n * n * 3];
    for j in 0..n {
        // png rows run top to bottom; the field stores j = 0 at the bottom.
        let dst = n - 1 - j;
        for i in 0..n {
            let t = ((field[j * n + i] - lo) / range).clamp(0.0, 1.0);
            let [r, g, b] = colormap(t);
            let base = (dst * n + i) * 3;
            px[base] = r;
            px[base + 1] = g;
            px[base + 2] = b;
        }
    }
    let mut bytes = Vec::new();
    {
        let mut enc = png::Encoder::new(&mut bytes, n as u32, n as u32);
        enc.set_color(ColorType::Rgb);
        enc.set_depth(png::BitDepth::Eight);
        let mut writer = enc.write_header().map_err(|e| format!("png header: {e}"))?;
        writer
            .write_image_data(&px)
            .map_err(|e| format!("png data: {e}"))?;
    }
    Ok(bytes)
}
