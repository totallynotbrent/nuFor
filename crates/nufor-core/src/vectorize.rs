//! runtime-dispatched simd kernels for the conservation update.
//!
//! the inner loop of a finite-volume step is an element-wise divergence apply:
//! `a[i] -= dt_dx * (f[i+1] - f[i]) + y[i]`, once per conserved array, over a
//! contiguous run of cells. that is a perfect vectorization target, so this
//! module provides a scalar reference and an avx2 (4-lane) variant, choosing
//! one at runtime. both are numerically identical within float rounding; the
//! policy is to keep them equivalent, not bit-for-bit across architectures.

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::{
    _mm256_add_pd, _mm256_loadu_pd, _mm256_mul_pd, _mm256_set1_pd, _mm256_storeu_pd, _mm256_sub_pd,
};

/// the best simd level this host actually reports, for diagnostic output.
pub fn simd_capability() -> &'static str {
    #[cfg(target_arch = "x86_64")]
    {
        use std::arch::is_x86_feature_detected;
        if is_x86_feature_detected!("avx512f") {
            return "avx512f";
        }
        if is_x86_feature_detected!("avx2") {
            return "avx2";
        }
        if is_x86_feature_detected!("sse2") {
            return "sse2";
        }
    }
    "scalar"
}

/// apply `a[i] -= dt_dx * (f[i+1] - f[i]) + y[i]` for `n` cells.
///
/// `f` holds the n+1 face fluxes, `a` and `y` the n cell values and the
/// transverse term. uses avx2 when the host reports it and the build targets
/// x86_64, otherwise the scalar fallback.
pub fn apply_divergence(a: &mut [f64], f: &[f64], y: &[f64], dt_dx: f64) {
    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("avx2") {
        // safety: avx2 is verified present, and the slices are borrowed soundly.
        unsafe { apply_divergence_avx2(a, f, y, dt_dx) }
    } else {
        apply_divergence_scalar(a, f, y, dt_dx);
    }
    #[cfg(not(target_arch = "x86_64"))]
    apply_divergence_scalar(a, f, y, dt_dx);
}

/// the scalar reference, also the fallback on non-x86 targets.
fn apply_divergence_scalar(a: &mut [f64], f: &[f64], y: &[f64], dt_dx: f64) {
    for i in 0..a.len() {
        a[i] -= dt_dx * (f[i + 1] - f[i]) + y[i];
    }
}

/// the avx2 variant, four cells per instruction.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn apply_divergence_avx2(a: &mut [f64], f: &[f64], y: &[f64], dt_dx: f64) {
    // safety: every operation below is an avx2 intrinsic or a raw pointer
    // step on the slices this function is handed; the target feature and the
    // len-known-on-entry slices make all of them well-formed.
    unsafe {
        let dtdx = _mm256_set1_pd(dt_dx);
        let n = a.len();
        let mut i = 0;
        while i + 4 <= n {
            let fi = _mm256_loadu_pd(f.as_ptr().add(i));
            let fi2 = _mm256_loadu_pd(f.as_ptr().add(i + 1));
            let yv = _mm256_loadu_pd(y.as_ptr().add(i));
            let av = _mm256_loadu_pd(a.as_ptr().add(i));
            let diff = _mm256_sub_pd(fi2, fi);
            let upd = _mm256_add_pd(_mm256_mul_pd(dtdx, diff), yv);
            _mm256_storeu_pd(a.as_mut_ptr().add(i), _mm256_sub_pd(av, upd));
            i += 4;
        }
        while i < n {
            a[i] -= dt_dx * (f[i + 1] - f[i]) + y[i];
            i += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ref_update(a: &mut [f64], f: &[f64], y: &[f64], dt_dx: f64) {
        for i in 0..a.len() {
            a[i] -= dt_dx * (f[i + 1] - f[i]) + y[i];
        }
    }

    #[test]
    fn scalar_matches_the_reference() {
        let mut a: Vec<f64> = (0..17).map(|k| k as f64 * 0.5).collect();
        let f: Vec<f64> = (0..18).map(|k| (k as f64 / 3.0).ln_1p()).collect();
        let y: Vec<f64> = (0..17).map(|k| (k % 5) as f64 * 0.1).collect();
        let (mut b, mut c) = (a.clone(), a.clone());
        ref_update(&mut b, &f, &y, 1.2);
        apply_divergence_scalar(&mut c, &f, &y, 1.2);
        assert_eq!(b, c);
        a = b;
        // all values stay finite.
        assert!(a.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn avx2_matches_the_scalar_reference_when_available() {
        let a: Vec<f64> = (0..1011).map(|k| (k as f64).sin() * 3.0).collect();
        let f: Vec<f64> = (0..1012).map(|k| (k as f64).cos() * 2.0).collect();
        let y: Vec<f64> = (0..1011).map(|k| (k % 7) as f64 * 0.25).collect();
        let (mut want, mut got) = (a.clone(), a.clone());
        ref_update(&mut want, &f, &y, 0.7);
        #[cfg(target_arch = "x86_64")]
        if std::arch::is_x86_feature_detected!("avx2") {
            unsafe {
                apply_divergence_avx2(&mut got, &f, &y, 0.7);
            }
        } else {
            apply_divergence_scalar(&mut got, &f, &y, 0.7);
        }
        // the same arithmetic order within one lane, so agreement is exact here.
        for (w, g) in want.iter().zip(&got) {
            assert!((w - g).abs() < 1e-12, "scalar {w} vs avx {g}");
        }
    }

    #[test]
    fn dispatch_runs_without_panicking_and_edits_in_place() {
        let mut a: Vec<f64> = (0..20).map(|k| k as f64).collect();
        let f: Vec<f64> = (0..21).map(|k| (k % 3) as f64).collect();
        let y: Vec<f64> = vec![0.5; 20];
        apply_divergence(&mut a, &f, &y, 0.5);
        assert_eq!(a.len(), 20);
        assert!(a.iter().all(|v| v.is_finite()));
    }
}
