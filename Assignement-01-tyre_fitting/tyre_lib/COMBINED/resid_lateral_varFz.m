function res = resid_lateral_varFz(P,FY,KAPPA,ALPHA,GAMMA,FZ,tyre_data)

  % ----------------------------------------------------------------------
  %% Compute the residuals - least squares approach - to fit the Fx curve 
  %  with Fz=Fz_nom, IA=0. Pacejka 1996 Magic Formula
  % ----------------------------------------------------------------------

  % Define MF coefficients

  %Fz0 = 200*4.44822; % Nominal load 200 lbf
  
  tmp_tyre_data = tyre_data;
  
  tmp_tyre_data.rVy2 = P(1);
  
  %dfz = (Z - Fz0)./Fz0 ;
  
  % Longitudinal Force (Pure Longitudinal Slip) Equations
  res = 0;
  for i=1:length(ALPHA)
    [~, Gyk, SVyk] =  MF96_FXFYCOMB_coeffs_eqns(KAPPA(i), ALPHA(i), GAMMA(i), FZ(i), tmp_tyre_data);
    fy  = Gyk * MF96_FY0(0, ALPHA(i), GAMMA(i), FZ(i), tmp_tyre_data) + SVyk;
    res = res+(fy-FY(i))^2;
      %res = res+(fx0/FX(i)-1)^2;
  end
  
  % Compute the residuals
  res = res/sum(FY.^2);

end

