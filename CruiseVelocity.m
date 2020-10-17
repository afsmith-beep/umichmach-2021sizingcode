function [v_cruise] = CruiseVelocity(thrust_to_weight, MTOW, v, air_density, wing_ref_area, CD_0, AR, span_wing)
%Function to solve for cruise velocity
%Definitely needs to be touched up

%Coefficients for a rough dynamic thrust curve
T_0 = thrust_to_weight*MTOW*0.7;
T_1 = -0.060;
T_2 = -0.015;

v_cruise = zeros(size(sensor));

for j = [1:size(sensor,1)]
    j;
    for i = [1:size(span_wing,2)]
        v_cruise(j,i) = double(max(vpasolve(v.^2*T_2+v*T_1+T_0(j,i) ==...
        0.5*air_density*v.^2*wing_ref_area(j,i).*(CD_0+1./(pi*e*AR(j,i))*...
        (2*MTOW(j,i)./(air_density*v^2*wing_ref_area(j,i))).^2)...
        +1/2*air_density*v.^2*wing_ref_area(j,i)*Cf*wing_ref_area(j,i)/...
        span_wing(j,i)+1/2*air_density*v^2*0.02508382*0.36,v)));

        if imag(v_cruise(j,i)) || real(v_cruise(j,i)) < 0
            v_cruise(j,i) = 0;
        end
        
    end
end
end

