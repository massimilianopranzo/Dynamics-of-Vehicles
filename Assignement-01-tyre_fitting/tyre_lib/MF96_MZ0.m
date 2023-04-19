% Self aligning moment MZ0
function [mz0] = MF96_MZ0(alpha, phi, Fz, FY, tyre_data)

 % precode

  [alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr] = MF96_MZ0_coeffs(alpha, phi, Fz, tyre_data);

 % main code

  mz0 = magic_formula_MZ(alpha, alpha__t, alpha__r, Bt, Br, Ct, Dt, Dr, Et, FY);
  
 end
