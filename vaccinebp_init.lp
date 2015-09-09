\ SCIP STATISTICS
\   Problem name     : vaccinebp
\   Variables        : 10 (4 binary, 0 integer, 0 implicit integer, 6 continuous)
\   Constraints      : 23
\   Obj. scale       : 1
\   Obj. offset      : 0
Minimize
 Obj: +1.76197314026085 ww_2 +2.23802685973915 ww_3 +1.76197314026085 ww_4 +2.23802685973915 ww_5 +1.76197314026085 ww_6
      +2.23802685973915 ww_7
Subject to
 sh2_0_0: +4 Z0 +4 ww_2 >= +3
 sh2_0_1: +4 Z1 +4 ww_2 >= +3
 sh2_1_0: +5 Z2 +5 ww_3 >= +4
 sh2_1_1: +5 Z3 +5 ww_3 >= +4
 sh3_0_0: +4 Z0 >= +3
 sh3_1_0: +5 Z1 +5 ww_5 >= +4
 sh3_2_0: +4 Z2 +4 ww_6 >= +3
 sh3_3_0: +5 Z3 +5 ww_7 >= +4
 z0: +1 Z0 >= +0
 z0: +1 Z0 -1 ww_4 >= +0
 z1: +1 Z1 >= +0
 z1: +1 Z1 >= +0
 z2: +1 Z2 >= +0
 z2: +1 Z2 >= +0
 z3: +1 Z3 >= +0
 z3: +1 Z3 >= +0
 AlphaCon: +0.194034334187542 Z0 +0.24645895087767 Z1 +0.24645895087767 Z2 +0.313047764057118 Z3 <= +0.25
 w2_0: +1 ww_2 = +1
 w2_1: +1 ww_3 = +1
 w3_0: +1 ww_4 = +1
 w3_1: +1 ww_5 = +1
 w3_2: +1 ww_6 = +1
 w3_3: +1 ww_7 = +1
Bounds
 0 <= Z0 <= 1
 0 <= Z1 <= 1
 0 <= Z2 <= 1
 0 <= Z3 <= 1
Binaries
 Z0 Z1 Z2 Z3
End
