% Magic formula - stiffness for self aligning moment
function [mf_res] = magic_formula_stiffness_MZ_tmp(alpha, alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr, Fy)

 % precode

  

 % main code

  t1 = (Br ^ 2);
  t2 = (alpha__r ^ 2);
  t5 = sqrt((t2 * t1 + 1));
  t8 = sin(alpha);
  t11 = Bt * alpha__t;
  t12 = atan(t11);
  t16 = atan(-t11 + (t11 - t12) * Et);
  t18 = cos(t16 * Ct);
  mf_res = -t8 / t5 * Dr + t8 * t18 * Fy * Dt;
  
 end
