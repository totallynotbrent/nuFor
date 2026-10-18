# the 1d isentropic-with-shock reference for the cdv nozzle
# (liou, aiaa 87-0355). solves the quasi-1d flow for a given
# chamber/exit pressure ratio pe/pt:
#   pe/pt > p2 (second critical): fully subsonic flow
#   p3 <= pe/pt <= p2: choked, normal shock inside the divergent part
#   pe/pt < p3 (third critical): choked, fully supersonic divergent part
import math

g = 1.4

def area(x):
    if x < 5.0:
        return 1.75 - 0.75 * math.cos((0.2 * x - 1.0) * math.pi)
    return 1.25 - 0.25 * math.cos((0.2 * x - 1.0) * math.pi)

def area_ratio_mach(M):
    # A/A* on an isentrope
    return (1.0/M)*((2/(g+1))*(1+(g-1)/2*M*M))**((g+1)/(2*(g-1)))

def mach_from_ar(ar, supersonic):
    lo, hi = (1.0, 15.0) if supersonic else (1e-9, 1.0)
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if supersonic:
            lo, hi = (mid, hi) if area_ratio_mach(mid) < ar else (lo, mid)
        else:
            lo, hi = (mid, hi) if area_ratio_mach(mid) > ar else (lo, mid)
    return 0.5*(lo+hi)

def p_over_pt(M):
    return (1 + (g-1)/2*M*M)**(-g/(g-1))

def m2_of(M1):
    return math.sqrt((1 + (g-1)/2*M1*M1)/(g*M1*M1 - (g-1)/2))

def p02_over_p01(M1):
    p1 = p_over_pt(M1)
    p2_p1 = 1 + 2*g/(g+1)*(M1*M1 - 1)
    return p1 * p2_p1 / p_over_pt(m2_of(M1))

def critical_ratios(a_exit_over_astar=1.5):
    me_sup = mach_from_ar(a_exit_over_astar, True)
    return p_over_pt(me_sup) * p02_over_p01(me_sup), p_over_pt(me_sup)

# nondim reservoir for the solver runs: pt=1, Tt=3, R=1 -> rho_t=1/3
RHO_T = 1.0/3.0
A_T = math.sqrt(g*3.0)

def flux_at(M, A):
    fac = 1 + (g-1)/2*M*M
    rho = RHO_T / fac**(1/(g-1))
    a_ = A_T / math.sqrt(fac)
    return A * rho * M * a_

def subsonic_M_at(x, mdot):
    lo, hi = 1e-9, 0.999999
    for _ in range(300):
        M = 0.5*(lo+hi)
        if flux_at(M, area(x)) < mdot: lo = M
        else: hi = M
    return 0.5*(lo+hi)

def solve_1d(pe_pt, a_exit_over_astar=1.5):
    p_second, p_third = critical_ratios(a_exit_over_astar)
    # first critical: the exit pressure of the all-subsonic isentrope
    # whose throat just reaches M=1. above it the nozzle never chokes.
    pe_first = p_over_pt(mach_from_ar(a_exit_over_astar, False))
    xs = [i*10.0/400 for i in range(401)]
    out = []
    if pe_pt > pe_first:
        # fully subsonic: bisect the flow rate so the exit pressure
        # matches pe on the subsonic isentrope.
        lo, hi = 1e-8, flux_at(0.999999, area(5.0))
        for _ in range(200):
            mid = 0.5*(lo+hi)
            if p_over_pt(subsonic_M_at(10.0, mid)) < pe_pt: hi = mid
            else: lo = mid
        mdot = 0.5*(lo+hi)
        for x in xs:
            M = subsonic_M_at(x, mdot)
            out.append((x, M, p_over_pt(M)))
        kind = f"subsonic (mdot={mdot:.4f})"
    elif pe_pt >= p_third:
        # choked, normal shock inside the divergent section: bisect
        # the shock station so the exit pressure matches pe.
        def exit_p_for_shock_at(xs_shock):
            M1 = mach_from_ar(area(xs_shock), True)
            a2star = area(xs_shock) / area_ratio_mach(m2_of(M1))
            M_e = mach_from_ar(area(10.0)/a2star, False)
            return p02_over_p01(M1) * p_over_pt(M_e)
        lo, hi = 5.0001, 10.0
        for _ in range(200):
            mid = 0.5*(lo+hi)
            if exit_p_for_shock_at(mid) > pe_pt: lo = mid
            else: hi = mid
        xs_shock = 0.5*(lo+hi)
        M1s = mach_from_ar(area(xs_shock), True)
        a2star = area(xs_shock) / area_ratio_mach(m2_of(M1s))
        for x in xs:
            A = area(x)
            if x <= xs_shock:
                M = mach_from_ar(A, False) if x < 5.0 else mach_from_ar(A, True)
                out.append((x, M, p_over_pt(M)))
            else:
                M = mach_from_ar(A/a2star, False)
                out.append((x, M, p02_over_p01(M1s)*p_over_pt(M)))
        kind = f"shocked @ x={xs_shock:.3f}"
    else:
        for x in xs:
            A = area(x)
            M = mach_from_ar(A, False) if x < 5.0 else mach_from_ar(A, True)
            out.append((x, M, p_over_pt(M)))
        kind = "supersonic"
    return out, kind, (p_second, p_third)

if __name__ == "__main__":
    for pe in [0.89, 0.75, 0.16]:
        sol, kind, crit = solve_1d(pe)
        _, _, crit2 = solve_1d(pe)
        p_second, p_third = crit2
        pe_first = p_over_pt(mach_from_ar(1.5, False))
        print(f"pe/pt={pe}: 1st={pe_first:.4f} 2nd={p_second:.4f} 3rd={p_third:.4f} -> {kind}")
        with open(f"analytic_{str(pe).replace('.','')}.csv", "w") as f:
            f.write("x,M,p/pt\n")
            for x, M, p in sol:
                f.write(f"{x:.4f},{M:.5f},{p:.5f}\n")
    print("csv written")
