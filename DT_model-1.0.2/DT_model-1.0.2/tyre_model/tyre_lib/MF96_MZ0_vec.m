% Pure self aligning moment MZ0
% this function remap the scalar function to its vectorial form
function [mz0_vec] = MF96_MZ0_vec_tmp(alpha_vec, phi_vec, Fz_vec, tyre_data)
  
  mz0_vec = zeros(size(alpha_vec));
  for i = 1:length(alpha_vec)

    % [fx0_vec(i)] = MF96_FX0(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);

   	% precode
   	% [alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr] = MF96_MZ0_coeffs(alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);

   	% % main code
   	% [mz0_vec(i)] = magic_formula_MZ(alpha_vec(i), alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr);

    %%
    [mz0_vec(i)] = MF96_MZ0(alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);

  end
  
 end
