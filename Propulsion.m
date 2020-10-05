function [reqPropWeigh,  weightPropulsion] = Propulsion (MTOW, weightWings, weightFuselage, RegConst, AR, weightWings, thrustToWeight)
    
    for k = 0:100
        thrust = thrust_to_weight*MTOW; %thrust (N) 
        reqPropWeight = (thrust*0.224809 - RegConst(1))/RegConst(2)*4.44822; %Question: Where does 0.22.... come from?  Why is RegConst(1) subtracted?
        err = sum(sum(abs(reqPropWeight - weightPropulsion))); %sum of absolute error
        fprintf('err %f\n',err);

        if err < 1e-8
            fprintf('thrust converged after %d\n',k);
            break
        end
        
       weightPropulsion = weightPropulsion + 0.1*((thrust*0.224809 -RegConst(1))/RegConst(2)*4.44822 - weightPropulsion) ;
    end