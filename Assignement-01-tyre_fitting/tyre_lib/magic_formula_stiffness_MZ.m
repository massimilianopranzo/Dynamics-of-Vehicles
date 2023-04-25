% Magic formula - stiffness for self aligning moment
function [mf_res] = magic_formula_stiffness_MZ(alpha, alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr, B, C, D, E, SV, Fz, tyre_data)

 % precode

  [Fy] = MF96_FY0(0, alpha, 0, Fz, tyre_data);

 % main code

  t1 = (Br ^ 2);
  t2 = (alpha__r ^ 2);
  t5 = sqrt((t2 * t1 + 1));
  t8 = sin(alpha);
  t11 = B ^ 2;
  t12 = alpha ^ 2;
  t20 = B * alpha;
  t21 = atan(t20);
  t24 = -t20 + (t20 - t21) * E;
  t25 = t24 ^ 2;
  t30 = atan(t24);
  t31 = t30 * C;
  t32 = cos(t31);
  t34 = Bt * alpha__t;
  t35 = atan(t34);
  t39 = atan(-t34 + (t34 - t35) * Et);
  t41 = cos(t39 * Ct);
  t42 = cos(alpha);
  t46 = sin(t31);
  mf_res = -t8 / t5 * Dr + t42 * t41 * Dt * t32 / (t25 + 0.1e1) * (-B + (B - 0.1e1 / (t12 * t11 + 0.1e1) * B) * E) * D * C + t8 * t41 * Dt * (-t46 * D + SV);
  
 end
