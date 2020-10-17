function [sensorLength,sensorWeight,sensorContainerWeight] = Payload(n,max_sensor_length, min_sensor_length,max_sensor_mass,min_sensor_mass,min_container_mass,max_container_mass)
% Calculates various weight and size factors directly related to required
% payload. Outputted weights are in Newtons. Outputted Lengths are in
% meters
%   2020-2021 competition requires that sensors be carried inside within
%   cargo containers. This function will evalute a range of length and mass
%   combinations for those sensors

%% =========================== Unit Conversions ============================= %%

in2m = 0.0254; %inches to meters
ft2m = in2m * 12; %feet to meters
kg2N = 9.81; % Convert kg to newtons
oz2kg = 0.0283495; % Convert ounces to kg

%% =========================== Payload Function ============================= %%
sensorLength = linspace(min_sensor_length,max_sensor_length,n); % Sensor length test case vector in inches
sensorLength = sensorLength * in2m; % Sensor length converted to meters
sensorWeight = linspace(min_sensor_mass,max_sensor_mass,n); % Sensor weight test case vector in ounces
sensorWeight = sensorWeight * oz2kg; % Sensor mass in kg
sensorWeight = sensorWeight * kg2N; % Convert sensor mass to Newtons
sensorContainerWeight = linspace(min_container_mass,max_container_mass,n) * oz2kg; % Sensor container mass in kg
sensorContainerWeight = sensorContainerWeight * kg2N; % Sensor container weight in newtons

end

