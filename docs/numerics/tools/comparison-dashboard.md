---
title: Comparison dashboard
---

The comparison dashboard overlays several probe curves on one chart so you can
compare how a scalar changes along a line across mesh resolutions and between
fields. It is the second half of the field-visualization pair that began with
line probes, and it is the "comparison between runs" view of the web UI.

## What it compares

A comparison is built from two independent axes:

- **fields** — density, Mach number, or pressure along the cut;
- **mesh resolutions** — the probe is re-solved and re-sampled at each grid
  size.

Every combination is a series, so asking for `rho` and `mach` at `32, 64, 128`
cells yields six overlaid curves. The same line preset and sample count is
used across all of them, which is what makes the comparison fair.

## The API

`GET /api/compare` takes the fields, mesh resolutions, line preset, and sample
count, and returns a labelled list of series:

```
/api/compare?fields=rho,mach&n=32,64&line=H&samples=12
{
  "line": "H",
  "series": [
    {"label":"rho · n=32","field":"rho","n":32,"points":[{"s":0.0,"v":1.0}, ...]},
    {"label":"rho · n=64","field":"rho","n":64,"points":[ ... ]},
    ...
  ]
}
```

Line presets are shared with the single-probe endpoint: `H` (horizontal
mid-plane), `V` (vertical mid-plane), or `D` (diagonal).

## The web panel

The Compare tab lists field and resolution choices, a line preset, and a
sample count. Draw fetches `/api/compare` and plots every series on one shared
SVG, each in its own color, with the field and resolution labelled in a legend
below. Because all curves share the same axes, refinement trends and the
relative magnitudes of different scalars are visible at a glance.

## Verification

- A server test drives `/api/compare?fields=rho,mach&n=32,64&line=H&samples=12`
  and checks the response is JSON, contains both fields and both resolutions,
  and carries exactly four labelled series with no invalid values.
- Live `curl` shows four well-formed series on the diagonal cut, each with the
  requested sample count and a physically sensible profile.

The dashboard and the probe share their interpolation core, so the existing
probe tests cover the numerics underneath the comparison.

## See also

- [[numerics/tools/line-probes||Line probes]]
- [[numerics/tools/web-ui||Web UI skeleton]]
- [[numerics/2d/2d-euler||2D Euler solver]]
- [[verification/regression-suite|Regression suite]]