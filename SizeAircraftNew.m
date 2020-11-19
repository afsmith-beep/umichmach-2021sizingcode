function [wing_ref_area, AR, thrust, MTOW, Cl_takeoff, weight_propulsion] = SizeAircraftNew(span_wing, wing_ref_area, num_wings, dens_lin_wing, weight_fuselage, weight_propulsion, sensorWeight, sensorContainerWeight, thrust_to_weight, RegConst, airfoil_Cl_max, delta_Cl, air_density, Takeoff_velocity, sensor)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
for i = 0:1000

    % calculate the wingspan from the aspect ratio and wing reference area
    %span_wing = sqrt(AR .* wing_ref_area); %wingspan (m)
    AR = span_wing.^2./wing_ref_area; %wingspan (m)
    % wing weight in Newtons 
    weight_wings = num_wings*(span_wing*dens_lin_wing);
    for k = 0:100
        weight_empty = weight_fuselage + weight_wings + weight_propulsion; %empty weight (N)
        MTOW = weight_empty + sensor .* (sensorWeight + sensorContainerWeight); %Max Takeoff weight (N) (pass weigh 3 oz 0.085 kg)
        thrust = thrust_to_weight*MTOW; %required thrust (N) from MTOW
        req_prop_weight = ((thrust*0.224809 -RegConst(1))/RegConst(2))*4.44822;
        err = sum(sum(abs(req_prop_weight - weight_propulsion))); %sum of absolute error
        fprintf('err %f\n',err);

        if err < 1e-8
            fprintf('thrust converged after %f\n',k);
            break

        end
        
        weight_propulsion = weight_propulsion + ...
        0.1*((thrust*0.224809 -RegConst(1))/RegConst(2)*4.44822 - weight_propulsion) ;

    end
    
    Cl_stall = airfoil_Cl_max * AR ./ (AR + 2); % finite wing correction 
    Cl_takeoff = Cl_stall/(1.1^2)+delta_Cl; % equation from 481
    wing_area_req = 2*MTOW./(Cl_takeoff .* air_density * (Takeoff_velocity^2)); %required wing area


    err_wing_ref_area = sum(sum(abs(wing_ref_area - wing_area_req)));
    fprintf('err  %d  %f\n ', i, err_wing_ref_area);

    if err_wing_ref_area < 1e-8
        fprintf('wing_ref_area converged after %d\n',i);
        break

    end

    wing_ref_area = wing_ref_area + 0.1*(wing_area_req - wing_ref_area);  %m^2
end
end

