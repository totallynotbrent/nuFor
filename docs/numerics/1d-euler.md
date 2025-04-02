---
title: 1D Euler
description: The first solver milestone — a conservative 1D Euler solver for an ideal gas.
tags: [numerics]
---
# 1D Euler

The first serious numerical milestone is a conservative, finite-volume 1D
Euler solver for an ideal gas, evolving density, momentum, and total energy
and recovering physical shocks, contacts, and expansions.

This page tracks the milestone and the decisions behind it. It is the current
focus of the project.

## What it needs

- Conservative and primitive state and the conversion between them.
- Ideal-gas equation of state with pressure recovery and explicit
  physical-validity checks.
- Uniform 1D grid geometry, face-flux storage, first-order reconstruction,
  one baseline approximate flux (HLL, then HLLC), CFL time-step control, and
  explicit time advancement.
- Boundaries, residual computation, iteration logging, and restart writing.

## Verification ladder

Constant state → uniform advection/contact → Sod tube → Lax tube →
stationary shock → isentropic expansion.

A feature ships with a focused verification case and evidence; a successful
run or a screenshot is not verification.

## Related

- [[architecture|Architecture]]
- [[case-toml|Formats — case.toml (case definition)]]
- [[ffi-boundary|Research — FFI boundary]]