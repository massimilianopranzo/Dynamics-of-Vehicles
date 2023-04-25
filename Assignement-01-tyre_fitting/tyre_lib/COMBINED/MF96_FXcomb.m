% Pure longitudinal force FX0
function [fx] = MF96_FXcomb(kappa, alpha, phi, Fz, tyre_data)

  len = length(kappa);
  fx = zeros(len, 1);
   % precode
  for i=1:len
    [kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(kappa(i), alpha(i), phi(i), Fz(i), tyre_data);
    [Gxa, ~, ~] =  MF96_FXFYCOMB_coeffs_eqns(kappa(i), alpha(i), phi(i), Fz(i), tyre_data);
   % main code
  
    fx0 = magic_formula(kappa__x, Bx, Cx, Dx, Ex, SVx);
  
    fx(i) = fx0*Gxa;
  end
end
