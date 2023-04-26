% Combined lateral force FY
function [fy] = MF96_FYcomb(kappa, alpha, phi, Fz, tyre_data)

  len = length(kappa);
  ones_vec = ones(len, 1);
  fy = zeros(len, 1);
   % precode
  for i=1:len
    [alpha__y, By, Cy, Dy, Ey, SVy] = MF96_FY0_coeffs(0*ones_vec(i), alpha(i), phi(i), Fz(i), tyre_data);
    [~, Gyk, SVyk] =  MF96_FXFYCOMB_coeffs_eqns(kappa(i), alpha(i), phi(i), Fz(i), tyre_data);
   % main code
  
    fy0 = magic_formula(alpha__y, By, Cy, Dy, Ey, SVy);
  
    fy(i) = fy0*Gyk + SVyk;
  end
end
