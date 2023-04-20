p_val = -100:1:100;
for i = 1: length(p_val)
  v_min(i) = resid_pure_Fx_varGamma(p_val(i),FX_vec, KAPPA_vec, GAMMA_vec, tyre_coeffs.FZ0, tyre_coeffs);
end

figure()
plot(p_val, v_min)