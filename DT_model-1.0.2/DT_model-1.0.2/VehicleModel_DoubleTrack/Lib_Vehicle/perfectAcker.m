function [delta__fr, delta__fl] = perfectAcker(delta,vehicle_data)

    % ----------------------------------------------------------------
    %% Function purpose: compute steering angles for front wheels, with 
    %%                   perfect Ackerman steering model
    % ----------------------------------------------------------------
    
    % Load vehicle data
    Lf = vehicle_data.vehicle.Lf;        
    Lr = vehicle_data.vehicle.Lr;                 
    Wf = vehicle_data.vehicle.Wf;       
    delta_f0 = vehicle_data.front_wheel.delta_f0; % toe angle [deg]


    % Perfect Ackermann steering law
    delta__fr = delta_f0*pi/180 + atan((2*Lr+2*Lf)*tan(delta)/(Wf*tan(delta)+2*Lr+2*Lf));
    delta__fl = -delta_f0*pi/180 - atan((2*Lr+2*Lf)*tan(delta)/(Wf*tan(delta)-2*Lr-2*Lf));

end

