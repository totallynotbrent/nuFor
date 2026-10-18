---
title: Nozzle cold flow
---

Rung 7a: the axisymmetric converging-diverging nozzle, the first
case on the ladder that is an engine component rather than a
reentry shape. The validation case is Liou's CDV nozzle (AIAA
87-0355, the NPARC/WIND verification archive version): a
cosine area law from chamber to exit, run at three exit
pressure ratios that exercise every regime of 1D nozzle theory.

## Geometry

The area law, in the archive's units:

$$
A(x) = \begin{cases}
1.75 - 0.75\cos\left((0.2x - 1)\pi\right), \quad 0 \le x < 5 \\
1.25 - 0.25\cos\left((0.2x - 1)\pi\right), \quad 5 \le x \le 10
\end{cases}
$$

so A(0) = 2.5, A(5) = 1.0 (the throat), A(10) = 1.5. The
axisymmetric radius is $r(x) = \sqrt{A(x)/\pi}$: 0.892 at the
inlet, 0.564 at the throat, 0.691 at the exit. The wall enters
nuFor as a polygon body mask (61 contour points closed over the
domain top), the south boundary is the symmetry axis, and the
march is the inviscid axisymmetric Euler lane on a 160x48 grid.

## Boundary conditions

Two new boundary conditions arrived with this rung.

The **reservoir inflow** (west) holds a stagnation state and
expands it to the local face pressure:

$$
p_{ghost} = p_{int}, \quad
T_{ghost} = T_t \left(\frac{p_{int}}{p_t}\right)^{(\gamma-1)/\gamma}, \quad
\rho_{ghost} = \frac{p_{int}}{R T_{ghost}}, \quad
u_{ghost} = \sqrt{2 c_p (T_t - T_{ghost})}
$$

so the inflow is not a fixed state: when the throat chokes and a
compression wave travels upstream, the ghost answers it with a
denser, slower state, the acoustic feedback a supersonic-inflow
ghost cannot give. The tangential velocity is zero; the ghost
never injects swirl. Getting that tangential component wrong
(one bug during development answered the tangential variable
with the expansion speed) rings the whole inlet plane with a
spurious radial momentum source, and the fix is what let the
mass flux settle.

The **subsonic back-pressure outflow** (east) copies the interior
but pins the ghost pressure to the case's back value, the
standard treatment for the pressure-ratio studies whose exit is
subsonic. The fully supersonic study uses the supersonic outflow.

The axisymmetric runner also gained restart chaining
(`[time] restart_in` / `restart_out`, same format as the planar
lane) so a long march can be split across invocations.

## The three regimes

With A_e/A_t = 1.5 and gamma 1.4 the 1D theory fixes the
critical ratios first: the first critical (all-subsonic flow
that just reaches M = 1 at the throat) sits at p_e/p_t = 0.8805,
the third critical (fully supersonic exit) at 0.1602, and the
shock-at-exit value at 0.1263.

- **p_e/p_t = 0.89** sits just above the first critical:
  subsonic throughout, the nozzle never chokes.
- **p_e/p_t = 0.75** is below it: the flow chokes and a normal
  shock stands in the divergent section where the post-shock
  branch reaches the back pressure. The 1D solution puts it at
  x = 7.56.
- **p_e/p_t = 0.16** is just below the third critical: fully
  supersonic downstream of the throat with the 1D exit Mach
  1.854.

The 1D reference (quasi-1D isentropic with the normal-shock
jump placed at the pressure-matching station) is the analytic
gate; the 2D axisymmetric solution differs from it by the radial
distribution the 1D theory folds into a single streamwise
state, so the centerline and the exit-plane profiles are the
comparison surfaces and the deviation at the wall is part of the
answer, not an error.

## Results

The runs are the deep-converged marches (t = 160 to 400, well
past the startup transient; the supersonic run's exit Mach is
identical between t = 160 and t = 400, so the steady state is
reached).

**p_e/p_t = 0.89 (subsonic).** The cleanest gate. Exit
centerline p/pt = 0.884 against the 0.890 back pressure (0.7%),
exit Mach 0.414 against the 1D 0.411 (0.7%), and the exit-core
mass flux within 3.6% of the 1D value. The centerline Mach
peaks at 0.51 at the throat and decelerates smoothly; the flow
never chokes, matching the regime call.

**p_e/p_t = 0.16 (fully supersonic).** The throat chokes
(centerline M = 1.017 at x = 5) and the expansion runs to the
exit with mass flux conserved to 0.06% across the inlet, throat
and exit planes, the tightest statement of steady state the rung
produces. The exit plane shows the expected axisymmetric
profile: M lowest on the axis (1.731, where the pressure is
highest), peaking mid-radius at 1.837 against the 1D 1.854, and
dropping in the wall-adjacent rim. Mass-averaged exit Mach 1.75,
about 5.5% below the 1D value. The shortfall is geometric: the
staircase wall passes an effective throat area about 4% under
the analytic one (0.957 against 1.0), and the whole solution
shifts with it.

**p_e/p_t = 0.75 (normal shock).** The flow chokes, accelerates
to M ~ 1.6 just past the throat, and a normal shock stands near
x = 5.75 to 6.2. The 1D prediction puts the shock at x = 7.56:
the 2D shock sits early, and the post-shock region carries
staircase-induced oscillation that the other two cases do not
show. The exit pressure still matches the back value (0.755
against 0.750), so the far field is right, but the shock station
and the post-shock wall pressure are the weakest part of this
rung, and the honest verdict is: captured, not yet resolved.
The cut-cell rung (13) is the ladder's answer to this.

## The staircase caveat

The wall is a staircase on a Cartesian mesh, so the wall
pressure carries the step-induced ripple the cut-cell and
unstructured rungs (11-13) are on the ladder to remove, the
effective throat area runs a few percent under the analytic
one, and the shock case's post-shock noise lives exactly where
the steps interact with a standing discontinuity. The
centerline and exit-plane gates are clean because the axis
region is far from the steps; the wall distribution is reported
as-is with that caveat attached.

![centerline Mach vs 1d](/research/nozzle-cdv-mach.png)

![wall pressure vs 1d](/research/nozzle-cdv-wallp.png)

## Gates summary

- Throat: the choked runs reach M = 1 at x = 5; the subsonic
  run stays below 1 everywhere.
- Exit state vs 1D: subsonic 0.7% on p and M; supersonic
  mass-averaged M within 5.5% with the peak tracking the 1D
  value; shock case exit pressure within 0.7% but the shock
  station early.
- Mass flux: 0.06% spread across planes on the supersonic run.
- Figures: `docs/research/nozzle-cdv-mach.png`,
  `nozzle-cdv-wallp.png`; gate numbers in
  `cases/nozzle-cdv/gates.json`.
