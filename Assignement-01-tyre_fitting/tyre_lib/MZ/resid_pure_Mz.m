function res = resid_pure_Mz(P,MZ,ALPHA,GAMMA,FZ,FY,tyre_data)

    % ----------------------------------------------------------------------
    %% Compute the residuals - least squares approach - to fit the Fx curve 
    %  with Fz=Fz_nom, IA=0. Pacejka 1996 Magic Formula
    % ----------------------------------------------------------------------

    % Define MF coefficients

    %Fz0 = 200*4.44822; % Nominal load 200 lbf
    
    tmp_tyre_data = tyre_data;
    
    tmp_tyre_data.qHz1  = P(1); 
    tmp_tyre_data.qBz1  = P(2);
    tmp_tyre_data.qCz1  = P(3);
    tmp_tyre_data.qDz1  = P(4);
    tmp_tyre_data.qEz1  = P(5);
    tmp_tyre_data.qEz4  = P(6);
    tmp_tyre_data.qBz9  = P(7);
    tmp_tyre_data.qDz6  = P(8);
    tmp_tyre_data.qBz10 = P(9);
    
   %dfz = (Z - Fz0)./Fz0 ;
    
    % Aligning moment Equations
    res = 0;
    for i=1:length(ALPHA)
       mz0  = MF96_MZ0(ALPHA(i), GAMMA, FZ, FY(i), tmp_tyre_data);
       res = res + (mz0-MZ(i))^2;
    end
    
    % Compute the residuals
    res = res / sum(MZ.^2);
    

end

