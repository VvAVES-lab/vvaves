"""M2 Water Below: bench prototype (Python, integer ledgers). Rain -> ground -> water table -> springs -> rivers -> sediment.
Profile of N columns along the downhill line (neighbours by index only, MC-9; every vector would take a z).
All stored amounts are INTEGER nanometres of water/solid per unit area, so every move is exact and ledgers balance to 0.
Material numbers are stand-ins: nature.chem has no real permeability/solubility yet (M5 finding), so they live in MATERIALS below,
to be moved into each material's .chem book."""
import random, math
NM = 10**9
# per material: porosity, K (m/s, hydraulic conductivity), erodibility (per hour, fraction of the capacity gap), solubility (m/yr of surface lowering in contact with water), insoluble if 0
MATERIALS = {
 'sandstone': dict(por=0.20, K=1e-4, ero=0.05, sol=0.0),
 'limestone': dict(por=0.05, K=1e-6, ero=0.01, sol=2e-4),
 'clay':      dict(por=0.45, K=1e-9, ero=0.10, sol=0.0),
}
YEAR = 365*86400
class Profile:
    def __init__(s, z, floor, mat, dx=150.0, geo_scale=1.0, evap_mm_day=2.0, manning=0.05, ct=0.6, settle=0.01):
        n = len(z); s.n, s.dx, s.geo, s.evap, s.man, s.ct, s.settle = n, dx, geo_scale, evap_mm_day, manning, ct, settle
        s.zb = [int(v*NM) for v in z]; s.zf = [int(v*NM) for v in floor]; s.mat = list(mat)
        s.void = [int(MATERIALS[m]['por']*(s.zb[i]-s.zf[i])) for i, m in enumerate(mat)]
        s.por0 = [MATERIALS[m]['por'] for m in mat]
        s.gw = [0]*n; s.sw = [0]*n; s.sed = [0]*n
        s.slope_out = max((z[-2]-z[-1])/dx, 1e-3)
        s.L = dict(rain=0, evap=0, out_sw=0, out_gw=0, out_sed=0, dissolved=0)   # ledger
        s.solid0 = sum(s.zb[i]-s.zf[i]-s.void[i] for i in range(n))
        s.spring = [0]*n; s.t = 0.0
    def thick(s, i): return s.zb[i]-s.zf[i]
    def por(s, i): return s.void[i]/max(1, s.thick(i))
    def K(s, i): return MATERIALS[s.mat[i]]['K']*(s.por(i)/s.por0[i])**3     # cubic law: widening a crack opens it fast
    def head(s, i): return (s.zf[i] + s.gw[i]/max(s.por(i), 1e-9))/NM         # m
    def slope(s, i):
        j = i+1
        return s.slope_out if j >= s.n else (s.zb[i]-s.zb[j])/NM/s.dx
    # ---- rules (one named function each) ----
    def rule_rain(s, rain_mm):
        for i in range(s.n):
            a = int(rain_mm[i]*1e6); s.sw[i] += a; s.L['rain'] += a
    def rule_evaporation(s, dt):
        e = int(s.evap*1e6*dt/86400)
        for i in range(s.n):
            m = min(s.sw[i], e); s.sw[i] -= m; s.L['evap'] += m
    def rule_infiltration(s, dt):
        for i in range(s.n):
            Ki = s.K(i) * 5.0                    # surface layer more open than bulk: x5
            m = min(s.sw[i], int(Ki*dt*NM), max(0, s.void[i]-s.gw[i])); s.sw[i] -= m; s.gw[i] += m
    def rule_darcy(s, dt):                        # q = K h_sat dH/dx, moves between neighbours
        n = s.n; out = [0]*n; mv = [0]*n          # mv[i] = depth moved i -> i+1 (nm of water per area)
        H = [s.head(i) for i in range(n)]; hs = [s.gw[i]/max(s.por(i),1e-9)/NM for i in range(n)]
        for i in range(n-1):
            Km = 2/(1/s.K(i)+1/s.K(i+1)); h = 0.5*(hs[i]+hs[i+1])
            q = Km*h*(H[i]-H[i+1])/s.dx            # m2/s
            mv[i] = int(q*dt/s.dx*NM)
        K = s.K(n-1); qo = K*hs[n-1]*hs[n-1]/s.dx   # seepage face at the outlet, head = aquifer floor
        mo = int(qo*dt/s.dx*NM)
        for i in range(n):                        # clamp: never move more than half of what is stored
            o = max(mv[i], 0) + max(-mv[i-1], 0) if i > 0 else max(mv[i], 0)
            o += mo if i == n-1 else 0
            lim = s.gw[i]//2
            if o > lim and o > 0:
                f = lim/o
                if mv[i] > 0: mv[i] = int(mv[i]*f)
                if i > 0 and mv[i-1] < 0: mv[i-1] = int(mv[i-1]*f)
                if i == n-1: mo = int(mo*f)
        for i in range(n-1): s.gw[i] -= mv[i]; s.gw[i+1] += mv[i]
        s.gw[n-1] -= mo; s.L['out_gw'] += mo
    def rule_springs(s):                           # table above the surface -> water comes out
        for i in range(s.n):
            ex = s.gw[i]-s.void[i]
            if ex > 0: s.gw[i] -= ex; s.sw[i] += ex; s.spring[i] += ex
    def rule_surface_flow(s, dt):                 # Manning: v = d^(2/3) S^(1/2) / n; water carries its sediment
        n = s.n; mw = [0]*n; ms = [0]*n
        for i in range(n):
            if s.sw[i] <= 0: continue
            S = s.slope(i)
            if S <= 0: continue
            d = s.sw[i]/NM; v = d**(2/3)*S**0.5/s.man
            f = min(1.0, v*dt/s.dx); mw[i] = int(s.sw[i]*f); ms[i] = s.sed[i]*mw[i]//s.sw[i]
        for i in range(n):
            if mw[i] == 0: continue
            s.sw[i] -= mw[i]; s.sed[i] -= ms[i]
            if i+1 < n: s.sw[i+1] += mw[i]; s.sed[i+1] += ms[i]
            else: s.L['out_sw'] += mw[i]; s.L['out_sed'] += ms[i]
    def rule_erode_deposit(s, dt):               # capacity ct*v*S*d (volume fraction x depth); erode below it, drop above it
        for i in range(s.n):
            d = s.sw[i]/NM; S = max(s.slope(i), 0.0)
            if d <= 0 or S <= 0: cap = 0
            else: cap = int(s.ct*(d**(2/3)*S**0.5/s.man)*S*d*NM)
            if s.sed[i] < cap and i < s.n-1:       # last column = fixed base level (a rock sill)
                e = int((cap-s.sed[i]) * min(1.0, MATERIALS[s.mat[i]]['ero']*s.geo*dt/3600))
                e = min(e, s.thick(i)//20)
                if e > 0:
                    r = s.void[i]*e//s.thick(i)                 # voids leave with the rock, solids go into the water
                    s.zb[i] -= e; s.void[i] -= r; s.sed[i] += e-r
                    if s.gw[i] > s.void[i]: ex = s.gw[i]-s.void[i]; s.gw[i] -= ex; s.sw[i] += ex
            elif s.sed[i] > cap:
                f = 1.0 if d <= 0 else min(1.0, s.settle*dt/d)
                m = int((s.sed[i]-cap)*f); s.sed[i] -= m; s.zb[i] += m
    def rule_dissolve(s, dt):                    # soluble rock in contact with water: surface lowers, rock opens up underground
        for i in range(s.n):
            sol = MATERIALS[s.mat[i]]['sol']
            if sol == 0: continue
            if s.sw[i] > 0 and i < s.n-1:
                d = min(int(sol*s.geo*dt/YEAR*NM), s.thick(i)//20)
                r = s.void[i]*d//s.thick(i); s.zb[i] -= d; s.void[i] -= r; s.L['dissolved'] += d-r
            if s.gw[i] > 0:
                d = min(int(sol*s.geo*dt/YEAR*NM*s.gw[i]/max(1,s.void[i])), max(0, int(0.5*s.thick(i))-s.void[i]))
                s.void[i] += d; s.L['dissolved'] += d
    def step(s, dt, rain_mm):
        s.rule_rain(rain_mm); s.rule_evaporation(dt); s.rule_infiltration(dt); s.rule_darcy(dt); s.rule_springs()
        s.rule_surface_flow(dt); s.rule_erode_deposit(dt); s.rule_dissolve(dt); s.t += dt
    # ---- ledgers: must be exactly 0 ----
    def water_error(s): return s.L['rain']-(s.L['evap']+s.L['out_sw']+s.L['out_gw']+sum(s.sw)+sum(s.gw))
    def solid_error(s):
        solid = sum(s.zb[i]-s.zf[i]-s.void[i] for i in range(s.n))
        return s.solid0-(solid+sum(s.sed)+s.L['out_sed']+s.L['dissolved'])

def storms(seed, n, dt, mm_year=1200.0, mean_gap_h=60.0, mean_len_h=3.0):
    """Stand-in for Engine.Weather rain output: seeded storm events, mm per step per column (uniform over the strip)."""
    rng = random.Random(seed); steps_h = dt/3600
    mean_rate = mm_year/(365*24)*(mean_gap_h+mean_len_h)/mean_len_h
    left = 0.0; gap = rng.expovariate(1/mean_gap_h); rate = 0
    while True:
        if left > 0: left -= steps_h; r = rate*steps_h
        else:
            gap -= steps_h; r = 0.0
            if gap <= 0: left = rng.expovariate(1/mean_len_h); rate = mean_rate*rng.uniform(0.5, 1.5); gap = rng.expovariate(1/mean_gap_h)
        yield [r]*n

def hillslope(n=40, dx=150.0, mat='sandstone', top=210.0, s_surface=0.03, thin=0.0075, flat_tail=8):
    """Valley profile: surface 3% slope, flat foot at the end; aquifer floor steeper so the aquifer thins downvalley."""
    z = []; cur = top
    for i in range(n):
        z.append(cur); cur -= s_surface*dx if i < n-flat_tail else 0.002*dx
    thick = [max(1.0, 40.0-thin*dx*i) for i in range(n)]
    return z, [z[i]-thick[i] for i in range(n)], [mat]*n

def run(years, dt=3600.0, seed=1, geo=1.0, mm_year=1200.0, mat='sandstone', n=40, **kw):
    z, f, m = hillslope(n=n, mat=mat); p = Profile(z, f, m, geo_scale=geo, **kw)
    g = storms(seed, n, dt, mm_year=mm_year)
    for _ in range(int(years*YEAR/dt)): p.step(dt, next(g))
    return p
