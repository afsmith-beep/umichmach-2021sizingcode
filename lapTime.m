function [laps_10min,t_3laps] = lapTime(v_cruise,lap_length)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

laps_10min = (v_cruise)./(lap_length/3.28) * 60 * 10; %laps in 10 minutes
t_3laps = 3*(lap_length/3.28)./v_cruise*0.8; %time for 3 laps

end

