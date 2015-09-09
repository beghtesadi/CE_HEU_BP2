\ SCIP STATISTICS
\   Problem name     : t_vaccinebp
\   Variables        : 72 (4 binary, 0 integer, 0 implicit integer, 68 continuous)
\   Constraints      : 23
\   Obj. scale       : 1
\   Obj. offset      : 0
Minimize
 Obj: +1.76197314026085 t_ww_2 +2.23802685973915 t_ww_3 +1.76197314026085 t_ww_4 +2.23802685973915 t_ww_5
      +1.76197314026085 t_ww_6 +2.23802685973915 t_ww_7 +55.0616606331516 Alpha2_0 +69.9383393668485 Alpha2_1
      +55.0616606331516 Alpha3_0 +41.9630036201091 Alpha3_1 +55.0616606331516 Alpha3_2 +41.9630036201091 Alpha3_3
      +55.0616606331516 Alpha2_0 +41.9630036201091 Alpha2_1 +19.5827350227176 Alpha3_3 +33.0369963798909 Alpha2_0
      +41.9630036201091 Alpha2_1 +33.0369963798909 Alpha2_0 +19.5827350227176 Alpha2_1 +33.0369963798909 Alpha2_0
      +19.5827350227176 Alpha2_1 +15.4172649772824 Alpha2_0 +15.4172649772824 Alpha2_0 +39.6443956558691 Alpha2_0
      +33.0369963798909 Alpha3_0 +19.5827350227176 Alpha3_1 +17.6197314026085 Alpha2_0 +17.6197314026085 Alpha2_0
      +1.76197314026085 Alpha2_0 +1.76197314026085 Alpha2_0 +1.76197314026085 Alpha2_0 +50.3556043441309 Alpha2_1
      +33.0369963798909 Alpha3_2 +22.3802685973915 Alpha2_1 +3.08345299545649 Alpha2_0 +3.08345299545649 Alpha2_0
      +3.08345299545649 Alpha2_0 +3.08345299545649 Alpha2_0 +2.23802685973915 Alpha2_1 +15.4172649772824 Alpha3_2
      +3.91654700454351 Alpha2_1 +3.91654700454351 Alpha2_1 +3.91654700454351 Alpha2_1 +3.5239462805217 Alpha2_0
      +3.5239462805217 Alpha2_0 +3.91654700454351 Alpha3_3 +3.91654700454351 Alpha3_3 +23.7866373935215 Alpha2_0
      +3.91654700454351 Alpha3_1 +3.91654700454351 Alpha3_1 +15.4172649772824 Alpha3_0 +23.7866373935215 Alpha2_0
      +2.20246642532606 Alpha2_0 +2.20246642532606 Alpha2_0 +39.6443956558691 Alpha3_2 +17.6197314026085 Alpha3_2
      +3.08345299545649 Alpha3_2 +22.0246642532606 Alpha3_2 +30.2133626064785 Alpha2_1 +17.6197314026085 Alpha3_0
      +3.08345299545649 Alpha3_0 +1.76197314026085 Alpha3_0
Subject to
 sh2_0_0: +4 t_Z0 +4 t_ww_2 +150 Alpha2_0 +100 Alpha2_0 +50 Alpha2_0 +50 Alpha2_0 +5 Alpha2_0 +10 Alpha2_0
  +10 Alpha2_0 +55 Alpha2_0 +6 Alpha2_0 >= +3
 sh2_0_1: +4 t_Z1 +4 t_ww_2 +150 Alpha2_0 +100 Alpha2_0 +50 Alpha2_0 +50 Alpha2_0 +5 Alpha2_0 +10 Alpha2_0
  +10 Alpha2_0 +55 Alpha2_0 +6 Alpha2_0 >= +3
 sh2_1_0: +5 t_Z2 +5 t_ww_3 +100 Alpha2_1 +50 Alpha2_1 +10 Alpha2_1 +55 Alpha2_1 >= +4
 sh2_1_1: +5 t_Z3 +5 t_ww_3 +100 Alpha2_1 +50 Alpha2_1 +10 Alpha2_1 +55 Alpha2_1 >= +4
 sh3_0_0: +4 t_Z0 +150 Alpha2_0 +150 Alpha3_0 +100 Alpha2_0 +50 Alpha2_0 +100 Alpha2_0 +100 Alpha3_0 +50 Alpha2_0
  +5 Alpha2_0 +10 Alpha2_0 +10 Alpha2_0 +50 Alpha3_0 +55 Alpha2_0 +6 Alpha2_0 +50 Alpha3_0 +10 Alpha3_0
  +5 Alpha3_0 >= +3
 sh3_1_0: +5 t_Z1 +5 t_ww_5 +150 Alpha2_0 +100 Alpha3_1 +100 Alpha2_0 +100 Alpha2_0 +50 Alpha2_0 +100 Alpha2_0
  +50 Alpha3_1 +50 Alpha2_0 +5 Alpha2_0 +5 Alpha2_0 +10 Alpha2_0 +10 Alpha2_0 +10 Alpha2_0 +10 Alpha3_1
  +55 Alpha2_0 +6 Alpha2_0 >= +4
 sh3_2_0: +4 t_Z2 +4 t_ww_6 +150 Alpha2_1 +150 Alpha3_2 +100 Alpha2_1 +50 Alpha2_1 +100 Alpha2_1 +100 Alpha3_2
  +50 Alpha2_1 +5 Alpha2_1 +50 Alpha3_2 +10 Alpha2_1 +100 Alpha3_2 +50 Alpha3_2 +10 Alpha3_2 +50 Alpha3_2
  >= +3
 sh3_3_0: +5 t_Z3 +5 t_ww_7 +150 Alpha2_1 +100 Alpha3_3 +100 Alpha2_1 +50 Alpha3_3 +50 Alpha2_1 +100 Alpha2_1
  +50 Alpha2_1 +5 Alpha2_1 +10 Alpha2_1 +10 Alpha3_3 >= +4
 z0: +1 t_Z0 -1 Alpha2_0 -1 Alpha2_0 -1 Alpha2_0 >= +0
 z0: +1 t_Z0 -1 t_ww_4 >= +0
 z1: +1 t_Z1 -1 Alpha2_0 -1 Alpha2_0 -1 Alpha2_0 >= +0
 z1: +1 t_Z1 >= +0
 z2: +1 t_Z2 >= +0
 z2: +1 t_Z2 >= +0
 z3: +1 t_Z3 >= +0
 z3: +1 t_Z3 >= +0
 AlphaCon: +0.194034334187542 t_Z0 +0.24645895087767 t_Z1 +0.24645895087767 t_Z2 +0.313047764057118 t_Z3
  <= +0.25
 w2_0: +1 t_ww_2 +1 Alpha2_0 +1 Alpha2_0 +1 Alpha2_0 +1 Alpha2_0 +1 Alpha2_0 +1 Alpha2_0 +1 Alpha2_0 +1 Alpha2_0
  +1 Alpha2_0 +1 Alpha2_0 +1 Alpha2_0 +1 Alpha2_0 +1 Alpha2_0 +1 Alpha2_0 +1 Alpha2_0 +1 Alpha2_0 +1 Alpha2_0
  +1 Alpha2_0 +1 Alpha2_0 +1 Alpha2_0 +1 Alpha2_0 +1 Alpha2_0 +1 Alpha2_0 +1 Alpha2_0 = +1
 w2_1: +1 t_ww_3 +1 Alpha2_1 +1 Alpha2_1 +1 Alpha2_1 +1 Alpha2_1 +1 Alpha2_1 +1 Alpha2_1 +1 Alpha2_1 +1 Alpha2_1
  +1 Alpha2_1 +1 Alpha2_1 +1 Alpha2_1 +1 Alpha2_1 +1 Alpha2_1 = +1
 w3_0: +1 t_ww_4 +1 Alpha3_0 +1 Alpha3_0 +1 Alpha3_0 +1 Alpha3_0 +1 Alpha3_0 +1 Alpha3_0 +1 Alpha3_0 = +1
 w3_1: +1 t_ww_5 +1 Alpha3_1 +1 Alpha3_1 +1 Alpha3_1 +1 Alpha3_1 +1 Alpha3_1 = +1
 w3_2: +1 t_ww_6 +1 Alpha3_2 +1 Alpha3_2 +1 Alpha3_2 +1 Alpha3_2 +1 Alpha3_2 +1 Alpha3_2 +1 Alpha3_2 +1 Alpha3_2
  = +1
 w3_3: +1 t_ww_7 +1 Alpha3_3 +1 Alpha3_3 +1 Alpha3_3 +1 Alpha3_3 +1 Alpha3_3 = +1
Bounds
 0 <= t_Z0 <= 0
 0 <= t_Z1 <= 1
 0 <= t_Z2 <= 1
 0 <= t_Z3 <= 1
 0 <= t_ww_2 <= 4.87589485242978
 0 <= t_ww_3 <= 0
 0 <= t_ww_4 <= 0
 0 <= t_ww_5 <= 3.70150422149982
 0 <= t_ww_6 <= 1.46276845572893
 0 <= t_ww_7 <= 9.21296799705299
 0 <= Alpha2_0 <= 0
 0 <= Alpha2_1 <= 0
 0 <= Alpha3_0 <= 0
 0 <= Alpha2_0 <= 0
 0 <= Alpha2_1 <= 0.230324199926325
 0 <= Alpha3_0 <= 0
 0 <= Alpha3_1 <= 0.231600379969016
 0 <= Alpha3_2 <= 0.0468085905833259
 0 <= Alpha3_3 <= 0.921296799705298
 0 <= Alpha2_0 <= 0
 0 <= Alpha2_1 <= 0.921296799705298
 0 <= Alpha2_0 <= 0
 0 <= Alpha2_1 <= 0.921296799705298
 0 <= Alpha2_0 <= 0
 0 <= Alpha2_0 <= 0
 0 <= Alpha2_0 <= 0
 0 <= Alpha2_0 <= 0
 0 <= Alpha2_0 <= 0
 0 <= Alpha3_0 <= 0
 0 <= Alpha2_0 <= 0
 0 <= Alpha2_0 <= 0
 0 <= Alpha2_0 <= 0
 0 <= Alpha2_1 <= 0
 0 <= Alpha2_1 <= 0
 0 <= Alpha2_0 <= 0
 0 <= Alpha2_0 <= 0
 0 <= Alpha2_0 <= 0
 0 <= Alpha2_0 <= 0
 0 <= Alpha2_1 <= 0
 0 <= Alpha2_0 <= 0
 0 <= Alpha2_0 <= 0
 0 <= Alpha3_0 <= 0
 0 <= Alpha2_1 <= 0
 0 <= Alpha3_0 <= 0
Binaries
 t_Z0 t_Z1 t_Z2 t_Z3
End
