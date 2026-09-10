---
title: Line probes
---

Line probes sample a 2D scalar field along a straight segment and chart the
value against distance along the line, turning a colormap into a quantitative
profile. They are the first of the field-visualization additions to the web
UI; a run comparison dashboard follows.

## The probe

`probe_line` takes a field array, the grid it lives on, two endpoints, and a
sample count. It walks the segment at `n` evenly spaced points and
bilinearly interpolates the field from the four surrounding cell centers at
each one, returning `(distance, value)` pairs. Points that fall outside the
domain are clamped to the boundary cell, so a probe that starts or ends past
the grid reads the edge value instead of failing.

Bilinear interpolation is exact for fields that are linear in `x` and `y`,
which the unit test checks directly: probing `f = x + 2y` on the interior
reproduces the analytic value to machine precision.

## The API

`GET /api/probe` computes a resolved blast wave, extracts the requested
scalar (rho, mach, or p), and samples it along the segment given by
`x0,y0,x1,y1` at `samples` points:

```
/api/probe?n=96&field=mach&x0=0.1&y0=0.2&x1=0.9&y1=0.8&samples=25
{
  "field": "mach",
  "samples": [{"s": 0.0, "v": 0.11568}, ...]
}
```

`s` is the distance from the start of the segment, so a horizontal cut
across the unit domain runs `s` from 0 to 1.

## The web panel

The Probe tab offers a field selector, a line preset (horizontal mid-plane,
vertical mid-plane, or diagonal), a sample count, and a Sample button. The
browser fetches `/api/probe`, auto-scales both axes, and draws the value
curve as an SVG polyline, with the min/max and sample count in the status
line. Changing the field or the line re-samples immediately.

## Verification

- Unit tests cover exactness on a linear field, monotonic behaviour of a
  row-major field along a mid-plane cut, and clamping at out-of-domain
  points.
- A server test drives `/api/probe?n=64&field=rho&samples=9` and checks the
  response is JSON with the expected content and a full set of finite,
  non-NaN samples.
- Live `curl` against the running server returns valid, parseable JSON with
  physically sensible density across a supernova blast profile.

This mirrors the existing pattern: a feature ships with its API route, its
panel, and its test in the same change.

## See also

- [[numerics/tools/web-ui||Web UI skeleton]]
- [[numerics/2d/2d||2D foundations]]
- [[numerics/2d/2d-euler||2D Euler solver]]
- [[numerics/2d/oblique-shock||Oblique shock]]