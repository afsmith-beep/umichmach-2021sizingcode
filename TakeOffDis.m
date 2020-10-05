function [takeoff_dist] = TakeOffDis (AR, mu, thrust_to_weight, MTOW, air_density, Takeoff_velocity, wing_ref_area, Cl_takeoff, CD_0 )
    % Takeoff Distance: Raymer 487
    K = 1 ./ (pi * e * AR);
    %Kt = thrust_to_weight - mu
    Kt = thrust ./ MTOW - mu;
    Ka = air_density ./ 2 ./ (MTOW ./ wing_ref_area) .* (mu .* Cl_takeoff - CD_0 - K  .* Cl_takeoff .^2);

    takeoff_dist = (1 / 2 / g ./ Ka .* log((Kt + Ka .* Takeoff_velocity .^2) ./ Kt)) * 3.28084; % answer in ft converted from m

end