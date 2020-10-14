% ===========================================================
%                                         .:,:`           
%                                      `:,:               
%     ...       ...     `;;;          :,,     ;;;     ,;; 
%     ,,,:     ,,,,     ####        `,,:      ###     ;## 
%     ,,,,,   .,,,,     #####      `,,:       ###     ;## 
%     ,,,,,. `,,,,,    +##+##      ,,,        ###     ;## 
%     ,,,,,, ,,,,,,    ### ###    ,,,.        ###     ;## 
%     ,,, ,,,,, ,,,   .##` ###   ,,,,  #+     ########### 
%     ,,,  ,,,  ,,,   ###  `##;  *######+.,,,,########### 
%     ,,,  `,.  ,,,   ###   ###   ::,       .,,:.     :## 
%     ,,,   `   ,,,  ##########.  :,,         '`,,,,, ;## 
%     ,,,       ,,,  ###########   ,,,        ##@.,,,,, # 
%     ,,,       ,,, .##;     ###    ,,        ###  ,,,,,` 
%     ,,,       ,,: ###       ###    :,       ###  :,,,,, 
%                                   `:                 
% ===========================================================
%
% SUMMARY: Scoring analysis for the 2020-2021 DBF competition

%% Initialization
% close all
clear all
close all
clc

%% =========================== Unit Conversions ============================= %%

in2m = 0.0254; %inches to meters
ft2m = in2m * 12; %feet to meters
kg2N = 9.81; % Convert kg to newtons
oz2kg = 0.0283495; % Convert ounces to kg

%% =========================== Variables and Constants ===================================== %%
n = 30; % Number of test cases to run
g = 9.81; % Accleration due to gravity in m/s^2

% Wing Properties
span_wing = 54 * in2m; % Wingspan converted from inches to meters. 2021 competition states max 5ft wingspan. 54 inches accounts for 6 inch fuselage thickness
num_wings = 1; % Number of wings

% Fuselage Properties
weight_fuselage_initial = 23; % Empty fuselage weight guess in Newtons (Does not include wings, payloads, or propulsion system)
weight_fuselage = weight_fuselage_initial; % Sets up place holder for fuselage weight iterations
weight_propulsion = 18; % inital guess for Propulsion System weight in Newtons

% Desired and Approximated Aerodynamic Properties
thrust_to_weight = 0.6; % Desired thrust to weight ratio
Takeoff_velocity = 14; % Desired takeoff velocity in m/s
CD_0 = 0.06; % Zero-lift drag coefficient guess (CFD model approximation would be better)
e = 0.80; % Oswald efficiency factor (Need better method of approximating)
mu = 0.02;  % Takeoff distance constant
airfoil_Cl_max = 1.46; % Maximum lift coefficient for chosen BOE103
delta_Cl = 0.6*cosd(-10); % delta cl due to flaps: Raymer 279, 0.6 = Ratio of flapped area and total area
Cf = 0.006; % Coefficient of skin friction (Need a serious evaluation over how feasible this is)

% Course Properties
lap_length = 4000 * ft2m; % Approximate course length in feet converted to meters
air_density = 1.12; % Air density in Tucson, AZ

% Symbolic Variables
v = sym('v','real');

%% ================ Payload ======================= %%
% To meet competition criteria, minimum Length/Diameter ratio for the
% sensor is 4 where the minimum diameter is 1 inch
n_sensors = 14; % Maximum number of sensors that our aircraft could feasible carry
max_sensor_length = 12; % Maximum sensor length value to evaluate in inches
min_sensor_length = 4; % Minimum sensor length value to evaluate in inches
min_sensor_mass = 8; % Minimum sensor mass in ounces (Would be better to get a linear density plot instead of a guess)
max_sensor_mass = 18; % Maximum sensor mass in ounces (Would be better to get a linear density plot instead)
min_container_mass = 6; % Minimum sensor container mass in ounces (Would be better to get a linear density plot instead)
max_container_mass = 12; % Maximum sensor container mass in ounces (Would be better to get a linear density model)

[sensor_length,sensor_weight,sensor_container_mass] = Payload(n, n_sensors,max_sensor_length, min_sensor_length,max_sensor_mass,min_sensor_mass,min_container_mass,max_container_mass);