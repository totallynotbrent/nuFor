//! equilibrium-air thermodynamic closure: the Srinivasan-Tannehill
//! curve fits of NASA RP-1181 (also NASA CR-181245 / ERI-Ames 86401),
//! subroutine TGAS1, transcribed from the CR-181245 FORTRAN listing
//! with a derivative-identity audit applied to every coefficient.
//!
//! independent variables of the fits:
//!   Y = log10(rho / 1.292)      rho in kg/m^3 (1.292 = sea-level density)
//!   Z = log10(e / 78408.4)       e in m^2/s^2 (78408.4 = R*T0, T0 = 273.15 K)
//! valid for 250 K to 15000 K per the report.
//!
//! GAMM (the effective gamma = p/(rho*e) + 1) is a piecewise fit over
//! three density regions (A: Y <= -4.5, B: -4.5 < Y <= -0.5, C: Y > -0.5),
//! each subdivided in Z:
//!   A: 0.65 | 1.5 | 2.2 | 3.05 | 3.4 (warn above 3.69)
//!   B: 0.65 | 1.5 | 2.22 | 2.95 | 3.4 (warn above 3.4)
//!   C: 0.65 | 1.7 | 2.35 | 2.9 (warn above 2.9)
//! Z <= 0.65 in each region is the calorically-perfect cold branch with
//! region-specific constants (A: 1.3965, B: 1.398, C: 1.3988).
//!
//! the fit form per block:
//!   GAMM = g1 + g2 + g3 + g4 + (g5 + g6 + g7 + g8) / DENO
//! with g1 = a0 + a1*Y, g2 = (b0 + b1*Y)*Z, g3 = (c0 + c1*Z + c2*Y)*Y*Y,
//! g4 = (d0 + d1*Y + d2*Z)*Z*Z, and DENO a Grabau sigmoid:
//!   region A col 1: DENO = 1 - exp(s), s = g9 as printed (already the exp arg form)
//!   all other blocks: DENO = 1 + exp(s), s = g9, clamped to +/- 30
//! the g5..g8 polynomial-only blocks (A5, C1) use GAMM = g1+g2+g3+g4.
//!
//! partials for the sound speed follow the listing:
//!   GAMMR = dGAMM/dY (log space), GAMME = dGAMM/dZ (log space)
//!   a^2 = e * ((GAMM-1)*(GAMM + GAMME/ln10) + GAMMR/ln10)
//! the listing derives a^2 via d(p)/d(rho) at fixed entropy; the fits
//! provide the log-space partials and the 1/ln10 converts them.
//!
//! transcription truth policy: coefficients came from the typed FORTRAN
//! listing, each audited against the listing's own derivative blocks
//! (GASnR/GASnE algebraic identities), cross-checked with the printed
//! Table A1/A2 headers for the Z breakpoints. the one OCR-corrupted
//! threshold (label 170) was resolved to 2.22 from the table headers.
//!
//! validated gates (see tests):
//! - cold air at 1 atm: p within 0.3% of 101325 Pa, a ~ 340 m/s
//! - M22 normal shock: rho jump ~16x > perfect-gas cap of 6x

use crate::Error;

/// ln(10), the log-space conversion constant from the listing.
const LN10: f64 = std::f64::consts::LN_10;
/// RHO0 of the fits: sea-level density in kg/m^3.
pub(crate) const RHO0: f64 = 1.292;
/// E0 of the fits: R*T0 in m^2/s^2 with T0 = 273.15 K.
pub(crate) const E0: f64 = 78408.4;
/// cold-branch GAMM constants per density region (A, B, C).
const COLD_GAMM: [f64; 3] = [1.3965, 1.398, 1.3988];
/// Z breakpoints per density region, inner-first.
const Z_BREAKS: [[f64; 5]; 3] = [
    [0.65, 1.5, 2.2, 3.05, 3.4],  // A
    [0.65, 1.5, 2.22, 2.95, 3.4], // B
    [0.65, 1.7, 2.35, 2.9, 0.0],  // C
];

/// one block of the GAMM fit: the polynomial parts and the sigmoid.
struct Block {
    /// g1 = a + b*Y.
    g1: (f64, f64),
    /// g2 = (a + b*Y)*Z.
    g2: (f64, f64),
    /// g3 = (a + b*Z + c*Y)*Y*Y.
    g3: (f64, f64, f64),
    /// g4 = (a + b*Y + c*Z)*Z*Z.
    g4: (f64, f64, f64),
    /// grabau numerator present (g5..g8).
    grabau: Option<Grabau>,
}

/// the g5..g8 numerator and its transition function.
struct Grabau {
    /// g5 = a + b*Y.
    g5: (f64, f64),
    /// g6 = (a + b*Y)*Z.
    g6: (f64, f64),
    /// g7 = (a + b*Z + c*Y)*Y*Y.
    g7: (f64, f64, f64),
    /// g8 = (a + b*Y + c*Z)*Z*Z.
    g8: (f64, f64, f64),
    /// sigmoid exponent s = a + b*Y + c*Z + d*Y*Z.
    s: (f64, f64, f64, f64),
    /// region A col 1 uses DENO = 1 - exp(s); others 1 + exp(s).
    minus: bool,
}

/// the twelve fit blocks in dispatch order:
/// A1..A5, B1..B4, C1..C3.
static BLOCKS: [Block; 12] = [
    Block {
        g1: (1.52792, -0.0126953),
        g2: (-0.613514, -0.0508262),
        g3: (-0.00549384, 4.7512e-05, -0.000318468),
        g4: (0.631835, 0.0334012, -0.219921),
        grabau: Some(Grabau {
            g5: (-49.6286, -11.7932),
            g6: (69.1028, 44.0405),
            g7: (5.09249, -1.40326, 0.208988),
            g8: (13.7308, -17.8726, -18.6943),
            s: (24.60452, -2.0, -20.93022, 0.0),
            minus: true,
        }),
    },
    Block {
        g1: (-17.0333, -0.508545),
        g2: (24.6299, 0.445617),
        g3: (-0.00895298, 0.00229618, -0.000289186),
        g4: (-11.0204, -0.0989727, 1.62903),
        grabau: Some(Grabau {
            g5: (18.6797, 0.519662),
            g6: (-24.1338, -0.434837),
            g7: (0.00916089, -0.00152082, 0.000346482),
            g8: (10.2035, 0.0970762, -1.3946),
            s: (-142.762, -1.647088, 76.60312, 0.8259346),
            minus: false,
        }),
    },
    Block {
        g1: (2.24374, 0.103073),
        g2: (-0.532238, -0.0559852),
        g3: (0.00356484, -0.000101359, 0.000159127),
        g4: (-0.0480156, 0.0106794, 0.0366035),
        grabau: Some(Grabau {
            g5: (-5.70378, -0.310056),
            g6: (5.01094, 0.180411),
            g7: (-0.00949361, 0.00194839, -0.000224908),
            g8: (-1.40331, -0.0279718, 0.120278),
            s: (113.9755, -4.985467, -42.23833, 2.009706),
            minus: false,
        }),
    },
    Block {
        g1: (-20.807, 0.40197),
        g2: (22.591, -0.2566),
        g3: (-0.00095833, 0.0023966, 0.00033671),
        g4: (-7.7174, 0.04606, 0.878),
        grabau: Some(Grabau {
            g5: (-217.37, -4.6927),
            g6: (181.01, 2.6621),
            g7: (-0.034759, 0.0064681, -0.00070391),
            g8: (-50.019, -0.38381, 4.5795),
            s: (454.4373, 12.50133, -137.6001, -3.641774),
            minus: false,
        }),
    },
    Block {
        g1: (-52.2951, -0.400011),
        g2: (45.6439, 0.224484),
        g3: (-0.00373775, 0.00243161, 0.000224755),
        g4: (-12.9756, -0.0279517, 1.22998),
        grabau: None,
    },
    Block {
        g1: (1.39123, -0.00408321),
        g2: (0.0142545, 0.0141769),
        g3: (0.000257225, 0.000652912, 8.46912e-05),
        g4: (0.062555, -0.00783637, -0.097872),
        grabau: Some(Grabau {
            g5: (5.80955, -0.182302),
            g6: (-9.62396, 0.179619),
            g7: (-0.0230518, 0.011872, -0.000335499),
            g8: (5.27047, -0.0365507, -0.919897),
            s: (14.2, 0.0, -10.0, 0.0),
            minus: false,
        }),
    },
    Block {
        g1: (-1.20784, -0.257909),
        g2: (5.02307, 0.287201),
        g3: (-0.00995577, 0.00523524, -0.000145574),
        g4: (-3.20619, -0.0750405, 0.651564),
        grabau: Some(Grabau {
            g5: (-6.62841, 0.0277112),
            g6: (7.30762, -0.076823),
            g7: (0.00719421, -0.00362463, 0.000162777),
            g8: (-2.33161, 0.0304767, 0.166856),
            s: (125.5324, 2.015335, -63.90747, -0.6515225),
            minus: false,
        }),
    },
    Block {
        g1: (-2.2646, -0.0782263),
        g2: (4.90497, 0.0718096),
        g3: (-0.00306443, 0.00174209, 2.84214e-05),
        g4: (-2.2475, -0.0131641, 0.333658),
        grabau: Some(Grabau {
            g5: (-14.7904, -0.176627),
            g6: (13.5036, 0.087728),
            g7: (-0.00213327, 0.000715487, 7.30928e-05),
            g8: (-3.95372, -0.00896151, 0.363229),
            s: (178.8542, 6.317894, -67.56741, -2.46006),
            minus: false,
        }),
    },
    Block {
        g1: (-16.6904, -0.258318),
        g2: (17.835, 0.154898),
        g3: (-0.00971263, 0.0039774, 9.043e-05),
        g4: (-5.94108, -0.0201335, 0.660432),
        grabau: Some(Grabau {
            g5: (85.469, 11.7554),
            g6: (-72.176, -7.15723),
            g7: (-0.041615, 0.0138147, 0.000545184),
            g8: (20.1758, 1.0899, -1.86438),
            s: (288.3262, 12.48536, -88.16985, -3.720309),
            minus: false,
        }),
    },
    Block {
        g1: (1.37062, 0.0129673),
        g2: (0.111418, -0.0326912),
        g3: (0.00106869, -0.00200286, 0.000238305),
        g4: (-0.106133, 0.0190251, 0.0030221),
        grabau: None,
    },
    Block {
        g1: (0.0343846, -0.233584),
        g2: (2.85574, 0.259787),
        g3: (-0.01089927, 0.00423659, 0.000385712),
        g4: (-1.94785, -0.0673865, 0.408518),
        grabau: Some(Grabau {
            g5: (-4.20569, 0.133139),
            g6: (4.51236, -0.166341),
            g7: (0.00167787, -0.00110022, 0.000306676),
            g8: (-1.35516, 0.0491716, 0.0752509),
            s: (175.7042, -2.163278, -88.33702, 1.897543),
            minus: false,
        }),
    },
    Block {
        g1: (-1.70633, -0.148403),
        g2: (4.23104, 0.13729),
        g3: (-0.00910934, 0.00385707, 0.000269026),
        g4: (-1.97292, -0.028183, 0.295882),
        grabau: Some(Grabau {
            g5: (34.158, -18.9972),
            g6: (-40.858, 13.0321),
            g7: (-0.801272, 0.275121, -0.000177969),
            g8: (16.0826, -2.23386, -2.08853),
            s: (256.8123, 173.7089, -90.5889, -58.38803),
            minus: false,
        }),
    },
];

/// GAMM and its log-space partials at (Y, Z).
pub(crate) fn gamm_and_partials(y: f64, z: f64) -> (f64, f64, f64) {
    // density region and z sub-block; block indices into BLOCKS follow
    // the dispatch order A1..A5 = 0..4, B1..B4 = 5..8, C1..C3 = 9..11.
    // each region's top block extends upward past its last break (the
    // listing warns there but still evaluates).
    let (region, blk) = if y <= -4.5 {
        if z <= Z_BREAKS[0][0] {
            (0usize, usize::MAX)
        } else if z <= Z_BREAKS[0][1] {
            (0, 0)
        } else if z <= Z_BREAKS[0][2] {
            (0, 1)
        } else if z <= Z_BREAKS[0][3] {
            (0, 2)
        } else if z <= Z_BREAKS[0][4] {
            (0, 3)
        } else {
            (0, 4)
        }
    } else if y <= -0.5 {
        if z <= Z_BREAKS[1][0] {
            (1usize, usize::MAX)
        } else if z <= Z_BREAKS[1][1] {
            (1, 5)
        } else if z <= Z_BREAKS[1][2] {
            (1, 6)
        } else if z <= Z_BREAKS[1][3] {
            (1, 7)
        } else {
            (1, 8)
        }
    } else if z <= Z_BREAKS[2][0] {
        (2usize, usize::MAX)
    } else if z <= Z_BREAKS[2][1] {
        (2, 9)
    } else if z <= Z_BREAKS[2][2] {
        (2, 10)
    } else {
        (2, 11)
    };
    if blk == usize::MAX {
        return (COLD_GAMM[region], 0.0, 0.0);
    }
    let b = &BLOCKS[blk];
    let g1 = b.g1.0 + b.g1.1 * y;
    let g2 = (b.g2.0 + b.g2.1 * y) * z;
    let g3 = (b.g3.0 + b.g3.1 * z + b.g3.2 * y) * y * y;
    let g4 = (b.g4.0 + b.g4.1 * y + b.g4.2 * z) * z * z;
    let mut gamm = g1 + g2 + g3 + g4;
    let mut gr;
    let mut ge;
    // analytic partials of the polynomial part
    // g1 = a + b*y; g2 = (a + b*y)*z; g3 = (a + b*z + c*y)*y^2; g4 = (a + b*y + c*z)*z^2
    gr = b.g1.1 + b.g2.1 * z + b.g3.2 * y * y + 2.0 * y * (b.g3.0 + b.g3.1 * z) + b.g4.1 * z * z;
    ge = b.g2.0
        + b.g2.1 * y
        + b.g3.1 * y * y
        + 2.0 * z * (b.g4.0 + b.g4.1 * y)
        + 3.0 * b.g4.2 * z * z;
    // (g4 = (d0 + d1*y + d2*z)*z^2; d/dZ = 2*z*(d0 + d1*y) + 3*d2*z^2)
    if let Some(gb) = &b.grabau {
        let g5 = gb.g5.0 + gb.g5.1 * y;
        let g6 = (gb.g6.0 + gb.g6.1 * y) * z;
        let g7 = (gb.g7.0 + gb.g7.1 * z + gb.g7.2 * y) * y * y;
        let g8 = (gb.g8.0 + gb.g8.1 * y + gb.g8.2 * z) * z * z;
        let num = g5 + g6 + g7 + g8;
        let numr = gb.g5.1
            + gb.g6.1 * z
            + gb.g7.2 * y * y
            + 2.0 * y * (gb.g7.0 + gb.g7.1 * z)
            + gb.g8.1 * z * z;
        let nume = gb.g6.0
            + gb.g6.1 * y
            + gb.g7.1 * y * y
            + 2.0 * z * (gb.g8.0 + gb.g8.1 * y)
            + 3.0 * gb.g8.2 * z * z;
        let s = gb.s.0 + gb.s.1 * y + gb.s.2 * z + gb.s.3 * y * z;
        let sr = gb.s.1 + gb.s.3 * z;
        let se = gb.s.2 + gb.s.3 * y;
        if gb.minus {
            // DENO = 1 - exp(s)
            let e_s = s.exp();
            let deno = 1.0 - e_s;
            gamm += num / deno;
            gr += (numr + num * e_s * sr / deno) / deno;
            ge += (nume + num * e_s * se / deno) / deno;
        } else {
            // DENO = 1 + exp(s), s clamped to +/- 30
            let s = s.clamp(-30.0, 30.0);
            let e_s = s.exp();
            let deno = 1.0 + e_s;
            gamm += num / deno;
            gr += (numr - num * e_s * sr / deno) / deno;
            ge += (nume - num * e_s * se / deno) / deno;
        }
    }
    (gamm, gr, ge)
}

/// pressure from (rho, internal energy per mass).
pub fn pressure_from_energy(rho: &[f64], e: &[f64]) -> Result<Vec<f64>, Error> {
    let n = rho.len();
    if n == 0 || e.len() != n {
        return Err(Error::InvalidArgs);
    }
    let mut p = Vec::with_capacity(n);
    for i in 0..n {
        let rho_ok = rho[i].is_finite() && rho[i] > 1e-12;
        let e_ok = e[i].is_finite() && e[i] > 1e-3;
        if !(rho_ok && e_ok) {
            // out-of-physical cell during a transient: clamp to a cold
            // low-density floor so the march can keep going; the next
            // step's fluxes will re-pressurize from real neighbors.
            let rho_c = if rho_ok { rho[i] } else { 1e-6 };
            let e_c = if e_ok { e[i].max(1.0e4) } else { 1.0e4 };
            p.push(0.4 * e_c * rho_c);
            continue;
        }
        let y = (rho[i] / RHO0).log10().clamp(-6.0, 1.0);
        let z = (e[i] / E0).log10().clamp(0.65, 3.4);
        let (g, _, _) = gamm_and_partials(y, z);
        let pv = (g - 1.0) * e[i] * rho[i];
        if pv.is_finite() && pv > 0.0 {
            p.push(pv);
        } else {
            p.push(0.4 * 1.0e4 * rho[i].max(1e-12));
        }
    }
    Ok(p)
}

/// sound speed from (rho, internal energy per mass).
pub fn sound_speed_from_energy(rho: &[f64], e: &[f64]) -> Result<Vec<f64>, Error> {
    let n = rho.len();
    if n == 0 || e.len() != n {
        return Err(Error::InvalidArgs);
    }
    let mut a = Vec::with_capacity(n);
    for i in 0..n {
        let rho_ok = rho[i].is_finite() && rho[i] > 1e-12;
        let e_ok = e[i].is_finite() && e[i] > 1e-3;
        if !(rho_ok && e_ok) {
            a.push(200.0);
            continue;
        }
        let y = (rho[i] / RHO0).log10().clamp(-6.0, 1.0);
        let z = (e[i] / E0).log10().clamp(0.65, 3.4);
        let (g, gr, ge) = gamm_and_partials(y, z);
        let asq = e[i] * ((g - 1.0) * (g + ge / LN10) + gr / LN10);
        a.push(if asq.is_finite() && asq > 0.0 {
            asq.sqrt()
        } else {
            200.0
        });
    }
    Ok(a)
}

/// internal energy per mass from (rho, pressure): inverts the fits by
/// bisection on Z (pressure is monotone in e at fixed rho over the
/// fits' window). used by case initialization, which stores (rho, u, p)
/// and needs e for the conserved state.
pub fn energy_from_pressure(rho: &[f64], p: &[f64]) -> Result<Vec<f64>, Error> {
    let n = rho.len();
    if n == 0 || p.len() != n {
        return Err(Error::InvalidArgs);
    }
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        if rho[i] <= 0.0 || p[i] <= 0.0 {
            return Err(Error::InvalidArgs);
        }
        // bracket e in [E0*10^-1, E0*10^4] (250 K to ~3e4 K equivalent)
        let (mut lo, mut hi) = (78408.4e-1, 78408.4e4);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            let p_mid = eqair_pressure_at(rho[i], mid);
            if p_mid < p[i] {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        out.push(0.5 * (lo + hi));
    }
    Ok(out)
}

/// the closure at a single point: pressure from (rho, e).
pub fn eqair_pressure_at(rho: f64, e: f64) -> f64 {
    if rho <= 0.0 || e <= 0.0 {
        return f64::NAN;
    }
    let y = (rho / RHO0).log10();
    let z = (e / E0).log10();
    let (g, _, _) = gamm_and_partials(y, z);
    (g - 1.0) * e * rho
}

/// the closure at a single point: sound speed from (rho, e).
pub fn eqair_sound_at(rho: f64, e: f64) -> f64 {
    if rho <= 0.0 || e <= 0.0 {
        return f64::NAN;
    }
    let y = (rho / RHO0).log10();
    let z = (e / E0).log10();
    let (g, gr, ge) = gamm_and_partials(y, z);
    let asq = e * ((g - 1.0) * (g + ge / LN10) + gr / LN10);
    if asq <= 0.0 {
        return f64::NAN;
    }
    asq.sqrt()
}

/// the tgas1 face-fit tables in flat form for the fortran bridge: 12 blocks
/// of polynomial + grabau coefficients, the region z-breaks, the cold-limit
/// gammas, and the (rho0, e0) normalization.
pub(crate) struct Tgas1Tables {
    pub g1: [(f64, f64); 12],
    pub g2: [(f64, f64); 12],
    pub g3: [(f64, f64, f64); 12],
    pub g4: [(f64, f64, f64); 12],
    pub g5: [(f64, f64); 12],
    pub g6: [(f64, f64); 12],
    pub g7: [(f64, f64, f64); 12],
    pub g8: [(f64, f64, f64); 12],
    pub s: [(f64, f64, f64, f64); 12],
    pub minus: [bool; 12],
    pub z_breaks: [[f64; 5]; 3],
    pub cold_gamm: [f64; 3],
    pub rho0: f64,
    pub e0: f64,
}

/// snapshot the tables (clone; called once per process at mpi init).
pub(crate) fn tgas1_tables() -> Tgas1Tables {
    let blocks = &BLOCKS;
    let mut minus = [false; 12];
    for (i, b) in blocks.iter().enumerate() {
        if let Some(g) = &b.grabau {
            minus[i] = g.minus;
        }
    }
    Tgas1Tables {
        g1: std::array::from_fn(|i| blocks[i].g1),
        g2: std::array::from_fn(|i| blocks[i].g2),
        g3: std::array::from_fn(|i| blocks[i].g3),
        g4: std::array::from_fn(|i| blocks[i].g4),
        g5: std::array::from_fn(|i| blocks[i].grabau.as_ref().map_or((0.0, 0.0), |g| g.g5)),
        g6: std::array::from_fn(|i| blocks[i].grabau.as_ref().map_or((0.0, 0.0), |g| g.g6)),
        g7: std::array::from_fn(|i| blocks[i].grabau.as_ref().map_or((0.0, 0.0, 0.0), |g| g.g7)),
        g8: std::array::from_fn(|i| blocks[i].grabau.as_ref().map_or((0.0, 0.0, 0.0), |g| g.g8)),
        s: std::array::from_fn(|i| {
            blocks[i]
                .grabau
                .as_ref()
                .map_or((0.0, 0.0, 0.0, 0.0), |g| g.s)
        }),
        minus,
        z_breaks: Z_BREAKS,
        cold_gamm: COLD_GAMM,
        rho0: RHO0,
        e0: E0,
    }
}
