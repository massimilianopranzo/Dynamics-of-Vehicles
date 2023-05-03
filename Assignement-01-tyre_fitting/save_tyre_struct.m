clear
close all
clc

tc = load("tyre_Hoosier_B1464.mat");
v=fieldnames(tc.tyre_coeffs);
names = string(length(v));
for i = 1 : length(v)
    names(i) = v{i};
    value(i) = tc.tyre_coeffs.(v{i});
end

write_latex_macro('min_values.tex', names, value, 'w')

fileID = fopen('optimization_table.tex','w');
formatSpec = '%s & & & & %6.3 %s \n';
for i=1:length(v)
  fprintf(fileID,'%s \t & \t & \t & \t & \t %6.3f \t %s \n', names(i), value(i), '\\');
end
fclose(fileID);