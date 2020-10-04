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

%% =========================== Varibles and Constants ===================================== %%
n = 30; % Number of test cases to run
g = 9.81; % Accleration due to gravity in m/s^2

% Wing Properties
span_wing = 54 * in2m; % Wingspan converted from inches to meters. 2021 competition states max 5ft wingspan. 54 inches accounts for 6 inch fuselage thickness
num_wings = 1; % Number of wings

% Fuselage Properties
weight_fuselage_intial = 23; % Empty fuselage weight guess in Newtons (Does not include wings, payloads, or propulsion system)
weight_fuselage = weight_fuselage_initial; % Sets up place holder for fuselage weight iterations

% Desired and Approximated Aerodynamic Properties
thrust_to_weight = 0.6; % Desired thrust to weight ratio
Takeoff_velocity = 14; % Desired takeoff velocity in m/s
CD_0 = 0.06; % Zero-lift drag coefficient guess (CFD model approximation would be better)
e = 0.80; % Oswald efficiency factor (Need better method of approximating)
airfoil_Cl_max = 1.46; % Maximum lift coefficient for chosen BOE103
delta_Cl = 0.6*cosd(-10); % delta cl due to flaps: Raymer 279, 0.6 = Ratio of flapped area and total area

% Course Properties
lap_length = 4000 * ft2m; % Approximate course length in feet converted to meters
air_density = 1.12; % Air density in Tucson, AZ

% Symbolic Variables
v = sym('v','real');