%% Perform some analysis on the data
extra_params = model_sim.extra_params;
states = model_sim.states;
% figure()
% sum_F = states.Fz_rr.Data + states.Fz_rl.Data + states.Fz_fr.Data + states.Fz_fl.Data;
% plot(states.Fz_rr.Time, sum_F)

% figure()
% w_f = 9.81*0.5*(vehicle_data.vehicle.m*vehicle_data.vehicle.Lr/(vehicle_data.vehicle.Lf+vehicle_data.vehicle.Lr));
% plot(states.Fz_fr.Time, states.Fz_fr.Data-w_f);
% hold on
% w_r = 9.81*0.5*(vehicle_data.vehicle.m*vehicle_data.vehicle.Lf/(vehicle_data.vehicle.Lf+vehicle_data.vehicle.Lr));
% plot(states.Fz_fr.Time,( states.Fz_rr.Data-w_r), '--');

figure()
plot(states.Fz_fr.Time, extra_params.Fy_fr.Data);
hold on 
plot(states.Fz_fr.Time, 1/(vehicle_data.vehicle.L)*(extra_params.Mz_fr.Data), '--');

