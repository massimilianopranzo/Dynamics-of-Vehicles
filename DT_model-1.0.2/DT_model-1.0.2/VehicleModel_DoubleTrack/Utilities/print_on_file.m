fileID = fopen('results_axle_char.txt','a');
fprintf(fileID,'%.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \n', ped_0(end), delta_D(end), mu_r(end), mu_f(end), alpha_f(end), alpha_r(end), Ay(end));
fclose(fileID);