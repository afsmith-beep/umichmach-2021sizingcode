function [score, M1, M2, M3] = Scoring(sensor, t_3laps, laps_10min, sensorLength, sensorWeight)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
M1 = zeros(size(sensor));
M1(t_3laps/60 < 5) = 1; % Mission 1 score 

M2 = (sensor./t_3laps); % mission 2 score
M2 = 1+M2./max(M2);

M3 = 2+(laps_10min .* sensorLength .* sensorWeight)./max(laps_10min .* sensorLength .* sensorWeight); % Mission 3 score
GM = 6./(sensor+5); % Ground Mission Score
score = (M1 + M2 + M3+GM);
end

