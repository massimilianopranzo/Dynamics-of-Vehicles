% Pure self align moment
% this function remap the sclar function to its vectorial form
function [Calfa_vec] = MF96_CorneringStiffness_MZ(alpha_vec, phi_vec, Fz_vec, FY_vec, tyre_data)


  % [alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr] = MF96_MZ0_coeffs(alpha, phi, Fz, tyre_data)
  
  fx0_vec = zeros(size(alpha_vec));
  for i = 1:length(alpha_vec)
   % precode
   [alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr] = MF96_MZ0_coeffs(alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);

   % main code
    Calfa_vec(i) = magic_formula_stiffness_MZ(alpha_vec(i), alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr, FY_vec(i));
  end
  
 end
