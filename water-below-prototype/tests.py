"""Done-tests for M2 Water Below. Run: python3 tests.py"""
import sys, os, time; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from water import *
res = []
def check(name, ok, detail): res.append(ok); print(('PASS ' if ok else 'FAIL ')+name+': '+detail)
t0 = time.time()
GEO = 200.0   # erosion/dissolving compression (like Living Rock's compressed healing time); water paths are NOT compressed
p = run(3, geo=GEO); z0, f0, _ = hillslope(); n = p.n
sat = [i for i in range(n) if p.gw[i] >= p.void[i]-1 or p.spring[i] > 0]
check('WB1 rain builds a water table', sum(p.gw) > 0 and p.gw[0] > 0 and p.head(0) > p.zf[0]/NM,
      f'stored groundwater {sum(p.gw)/1e9:.2f} m-column total, head under the ridge {p.head(0)-p.zf[0]/NM:.2f} m above the floor')
sp = [i for i in range(n) if p.spring[i] > 0]
upstream_dry = all(p.spring[i] == 0 for i in range(0, 10))
check('WB2 a spring appears where the table meets the surface', len(sp) > 0 and upstream_dry,
      f'springs at columns {sp[0] if sp else None}..{sp[-1] if sp else None} (x = {sp[0]*p.dx/1000 if sp else 0:.1f} km); none on the thick upper aquifer: {upstream_dry}')
dz = [(p.zb[i]-int(z0[i]*NM))/1e9 for i in range(n)]
chan = [i for i in range(10, n-8)]
mean_cut = -sum(dz[i] for i in chan)/len(chan)
check('WB3 the valley deepens over the years', mean_cut > 0 and min(dz) < -0.01 and abs(dz[0]) < abs(min(dz)),
      f'3 yr (erosion x{GEO:g}): mean lowering {mean_cut:.3f} m, deepest cut {min(dz):.3f} m at column {dz.index(min(dz))}, ridge top {dz[0]:+.4f} m')
tail = dz[n-8:n-1]
check('WB4 slow water drops sediment: a delta grows on the flat foot', max(tail) > 0 and sum(p.sed) >= 0,
      f'foot gained up to {max(tail):.3f} m; sediment exported {p.L["out_sed"]/1e9:.3f} m-column')
check('WB5 water and sediment ledgers balance to exactly 0', p.water_error() == 0 and p.solid_error() == 0,
      f'water error {p.water_error()} nm, solid error {p.solid_error()} nm (integers); rain in {p.L["rain"]/1e9:.2f}, evaporated {p.L["evap"]/1e9:.2f}, left via river {p.L["out_sw"]/1e9:.2f}, via seepage {p.L["out_gw"]/1e9:.2f}')
ls = run(3, geo=GEO, mat='limestone'); ss = p
k0 = [MATERIALS['limestone']['K']]*ls.n
gr=[ls.K(i)/MATERIALS['limestone']['K'] for i in range(ls.n)]; growth=max(gr); gi=gr.index(growth)
check('WB6 limestone dissolves into karst, sandstone does not', ls.L['dissolved'] > 0 and growth > 2 and growth > 2*gr[2] and ss.L['dissolved'] == 0 and ls.solid_error() == 0,
      f'limestone: dissolved {ls.L["dissolved"]/1e9:.3f} m-column, conductivity grew x{growth:.0f} at column {gi} (porosity {ls.por(gi):.3f}, from 0.05) vs x{gr[2]:.1f} under the ridge, so caves open where the water flows; sandstone dissolved {ss.L["dissolved"]}')
d = run(3, geo=GEO, mm_year=0.0)
check('WB7 no rain, no table, no spring, no erosion', sum(d.gw) == 0 and sum(d.spring) == 0 and d.zb == [int(v*NM) for v in z0], 'all zero')
a, b, c = run(1, seed=7), run(1, seed=7), run(1, seed=8)
check('WB8 same seed same history; other storms differ', a.zb == b.zb and a.gw == b.gw and a.gw != c.gw, 'identical rerun; seed 8 differs')
print(f'\n{sum(res)}/{len(res)} pass, {time.time()-t0:.0f} s total')
