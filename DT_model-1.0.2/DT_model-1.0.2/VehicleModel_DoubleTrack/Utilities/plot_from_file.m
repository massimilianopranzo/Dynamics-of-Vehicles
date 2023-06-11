data_read = readmatrix('results_axle_char.txt');
ped_0 = data_read(:,1);
delta_D = data_read(:,2);
mu_r = data_read(:,3);
mu_f = data_read(:,4);
alpfa_f = data_read(:,5);
alpfa_r = data_read(:,6);
Ay = data_read(:,7);

offset = 0.001;

%% Normalized axle characteristic
figure('Name','Norm axle char','NumberTitle','off'), clf
plot(alpfa_r, mu_r,  'o', 'LineWidth', 2)
hold on
plot(alpfa_f, mu_f,  'o', 'LineWidth', 2)
grid on
for i = 1:size(mu_r,1)
    text(mu_r(i), alpfa_r(i) + offset, num2str(Ay(i)))
end
legend('$\mu_{r}$', '$\mu_{f}$', 'location', 'southeast')