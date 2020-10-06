function [weightPropulsion] = Propulsion (MTOW, weightWings, weightFuselage, RegConst, AR, weightWings, thrustToWeight)
   % Estimated weight of the electronic systems 
   % Todo: calculate endurance based on electronic systems and battery
   % capacity
    for k = 0:100
        thrust = thrust_to_weight*MTOW; %thrust (N) 
        reqPropWeight = (thrust*0.224809 - RegConst(1))/RegConst(2)*4.44822; 
        err = sum(sum(abs(reqPropWeight - weightPropulsion))); %sum of absolute error
        fprintf('err %f\n',err);

        if err < 1e-8
            fprintf('thrust converged after %d\n',k);
            break
        end
        
       weightPropulsion = weightPropulsion + 0.1*((thrust*0.224809 -RegConst(1))/RegConst(2)*4.44822 - weightPropulsion) ;
    end