# M2 Water Below: bench prototype (2026-11-04)
Status: [tested] Python prototype in `prototype/` (water.py, tests.py, tests-2026-11-04.txt). Not yet in the repo. Next: port to `engine/Engine.Water` (C#) per MC-1..MC-9, R20/R21, once Cvpid replies in the thread and a Remote Control session can open.

## Shape
40-column profile along the downhill line (150 m columns), neighbours only by index (MC-9), nothing but SI amounts per column and a floor/bed elevation (R21). All stored amounts are integer nanometres of water or solid per area, so every move is exact and both ledgers balance to 0 (not just ~1e-13). Step 1 h. Rain comes in as mm per step per column (stand-in: seeded storms; the port reads Engine.Weather's rain accumulator).

## Rules (one named function each)
| Name | Equation | Note |
|---|---|---|
| rain / evaporation | sw += R; sw -= 2 mm/day | evaporation is a ledger sink |
| infiltration | m = min(sw, 5 K dt, free pore space) | surface layer x5 more open |
| darcy | q = K_mean h_sat (H_i - H_i+1)/dx, H = floor + stored/porosity | harmonic K, at most half a cell's store moves per step, seepage face at the outlet |
| springs | stored above pore space comes out as surface water | a spring is simply the table meeting the surface |
| surface_flow | v = d^(2/3) S^(1/2) / n (Manning); water carries its sediment | |
| erode_deposit | capacity = ct v S d; erode ero*(cap - load), drop settle*dt/d of any excess | last column is a fixed base level |
| dissolve | soluble rock under water lowers sol m/yr; groundwater in it widens voids | K = K0 (por/por0)^3, so caves open where water flows |

## Results
7 of 8 pass first run, 8 of 8 after fixing the karst test (below). Water error 0 nm, solid error 0 nm. 3 simulated years run in ~1.6 s.
## Honest notes
- Erosion and dissolving are compressed x200 (like Living Rock's healing time); water paths are not. Uncompressed, 3 years lowers the valley by micrometres. Real erosion is mm/kyr, so deep time needs the Day-lapse/deep-time machinery, not more steps.
- Karst check first failed: I tested the ridge, but caves open where water flows (the spring zone), not under the dry ridge. The test now checks growth is greatest where water flows (my correction of the test, not the rule).
- Material numbers (porosity, K, erodibility, solubility) are stand-ins; nature.chem has none of them (M5 finding). The port must put them into .chem books, or LC-4 cannot supply them.
- Sediment is one grain size and deposits instantly onto the bed; no bank collapse, no channel width, no 2D flow yet. Profile only: the column-world (CV-1) hook-up needs the 4/8-neighbour version plus lowest-neighbour routing.
- Not tested: rivers over years of drought/flood extremes, sensitivity of the constants (ct 0.6, settle 0.01, Manning 0.05).
## Prior art
Exists: Darcy groundwater, Manning flow, stream-power/transport-capacity erosion, karst dissolution feedback are textbook; games (From Dust, Minecraft) script rivers rather than derive springs. Ours: the pairing in tiny named rules with exact integer ledgers.
