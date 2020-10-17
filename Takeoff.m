function [wing_ref_area,takeoff_dist] = Takeoff(mu,e,AR,MTOW,thrust,air_density,wing_ref_area,Cl_takeoff,CD_0,g,Takeoff_velocity)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
%Takeoff Distance: Raymer 487

K = 1./(pi*e*AR);
Kt = thrust./MTOW-mu;
Ka = air_density./2./(MTOW./wing_ref_area).*(mu.*Cl_takeoff-CD_0-K.*Cl_takeoff.^2);

takeoff_dist = (1/2/g./Ka.*log((Kt+Ka.*Takeoff_velocity.^2)./Kt))*3.28084; % answer in ft converted from m
% MTOW = MTOW+( sensor * (0.141748 * g))
%change negative to a very large number and NaN values to a very small number
MTOW(MTOW<= 0) = 1e6;
MTOW(isnan(MTOW)==1) = 0.1;

%change negative to a very large number and NaN values to a very small number
wing_ref_area(wing_ref_area<= 0) = inf;
wing_ref_area(isnan(wing_ref_area)== 1) = 0.01;
end

