# Fay-Riddell stagnation heating

The rung-6 validation anchor: the closed form for stagnation-point
cold-wall heat flux on a blunt nose, implemented in
`nufor-core::fay_riddell` as the perfect-gas Lewis-number-one form.

## The form

$$q = 0.763\,\mathrm{Pr}^{-0.6} \sqrt{\rho_e \mu_e \beta}\;(h_{0e} - h_w)$$

with the Newtonian stagnation velocity gradient

$$\beta = \frac{1}{R_n}\sqrt{\frac{2(p_e - p_\infty)}{\rho_e}}$$

where the edge state is the post-normal-shock state at the nose and
$h_{0e}$ the total enthalpy behind that shock. Viscosity is
Sutherland's law for N2. Lewis number one means the
dissociation-enthalpy bracket of the full form is absent; that term
activates when the species rung lands, where the wall enthalpy picks
up the formation-enthalpy difference between edge and wall
composition.

## What is pinned

- The Rankine-Hugoniot tables at M=20, gamma=1.4: post-shock
  pressure ratio 465.5, density ratio 5.926, temperature ratio 78.55
  (all within 0.5 percent).
- The stagnation chain: normal-shock total-pressure ratio times the
  post-shock isentropic recovery, giving $p_0$, $T_0$, $h_{0e}$.
- Scaling: heating falls with the square root of nose radius, rises
  with the square root of edge density and viscosity, and vanishes
  with the wall-enthalpy difference. On the Mach-20 N2 cylinder case
  (0.89 Pa, 220 K, 1 m nose, 1000 K wall) the closed form reads
  $1.15 \times 10^5$ W/m², the expected order for that class.

## What is deliberately not claimed

The correlation constant in the Detra-Kemp-Riddell / Sutton-Graves
engineering family varies across secondary sources by orders of
magnitude depending on the assumed unit system, so no
cross-check gate against it is asserted here. The closed form's own
anchors are the shock tables and the scaling laws above; the
empirical comparison joins when a published table with stated units
is pinned.

## How the gate will close

Rung 6 runs a wall-resolved viscous sphere at entry conditions and
compares the solver's stagnation-line wall heat flux against this
closed form. The prerequisites are in the tree: the viscous lane
(Sutherland transport, isothermal wall), the body mask, and this
form. The remaining work is wall-normal resolution that holds a
boundary layer; the Blasius hold on the flat-plate viscous gate is
the resolution pattern.

The Euler-wall heat probe in the research pages is a grid artifact
by construction; nothing the solver prints about wall heating is
citable until this gate passes.
