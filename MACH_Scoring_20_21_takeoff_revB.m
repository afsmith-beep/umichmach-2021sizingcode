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
% clc

%% =========================== Unit Conversions ============================= %%

in2m = 0.0254; %inches to meters
ft2m = in2m * 12; %feet to meters

%% =========================== Varibles ===================================== %%

n = 30; % Number of test cases 
span_wing = linspace(54, 54, n)*in2m; % Generates vector of wingspans
sensor = linspace(1, 14, n); % Generates vector for number of sensors carried by aircraft
[span_wing, sensor] = meshgrid(span_wing, sensor); % Creates matrix relating cases for each wingspan/sensor configuration

g = 9.81; % m/s^2
num_wings = 1; %wings

weight_fuselage_intial = 23;  % fuse empty weight Newtons

weight_fuselage = weight_fuselage_intial; % Newtons 

%wing linear density - weight of the wing per inch
dens_lin_wing = (0.2875 * 0.0283495 * g / in2m);  % density N/m

thrust_to_weight = 0.6; % Desired aircraft thrust to weight ratio
Takeoff_velocity = 14; % Desired takeoff velocity in m/s
CD_0 = 0.06; % Zero lift drag coefficient 
e = 0.80; % oswald efficiency
lap_length = 4000; % Contest course length in ft
airfoil_Cl_max = 1.46; % Maximum lift coefficient (1.46)
delta_Cl = 0.6*cos(-10*pi/180); % delta cl due to flaps: Raymer 279, 0.6 = Ratio of flapped area and total area
% historical data
air_density = 1.12; % Tucson air density (kg/m^3)
%symbolic variables
v = sym('v','real');

wing_ref_area = 0.8; % inital guess for wing area (m^2)
weight_propulsion = 18; % inital guess for Propulsion System weight in Newtons
Cf=  0.006; % Skin friction coefficient estimate
RegConst = [-0.823411250229396;4.34939493516460];

sensorLength = linspace(4,12,n); % Sensor length test case vector
sensorLength = sensorLength * in2m; % Sensor length in meters
sensorWeight = linspace(8,18,n); % Sensor weight test case vector in ounces
sensorWeight = sensorWeight * 0.0283495; % Sensor mass in kg
sensorContainer = 6 * 0.0283495; % Sensor container mass in kg

%% ================ Data for Propulsion System Weight ======================= %%
%%%%%UPDATE%%%%%%%%
thrust_data = [0.5, 1, 2, 3, 4, 5, 15]*4.44822; % Thrust needed in Newtons
weight_data =[9, 12, 25.9, 37.5, 42, 55, 113]*0.28; % Weight of prop. sys. in N

prop_weight_coef = polyfit(thrust_data, weight_data, 1 );

%% ========================= Size Aircraft ================================= %%
% Iterative method to solve for wing_ref_area, weight_propulsion, MTOW, weight_empty
% Note: this method is not very sophisticated and prone to error 
weight_fuselage = 0.5*sensor+weight_fuselage; % This does NOT make sense
for i = 0:1000

    % calculate the wingspan from the aspect ratio and wing reference area
    %span_wing = sqrt(AR .* wing_ref_area); %wingspan (m)
    AR = span_wing.^2./wing_ref_area; %wingspan (m)
    % wing weight in Newtons 
    weight_wings = num_wings*(span_wing*dens_lin_wing);
    for k = 0:100
        weight_empty = weight_fuselage + weight_wings + weight_propulsion; %empty weight (N)
        MTOW = weight_empty + sensor .* ((sensorWeight + sensorContainer) .* g); %Max Takeoff weight (N) (pass weigh 3 oz 0.085 kg)
        % thrust = 70;
        thrust = thrust_to_weight*MTOW; %thrust (N) (weight_empty or MTOW????)
        req_prop_weight = (thrust*0.224809 -RegConst(1))/RegConst(2)*4.44822;
        err = sum(sum(abs(req_prop_weight - weight_propulsion))); %sum of absolute error
        fprintf('err %f\n',err);

        if err < 1e-8
            fprintf('thrust converged after %d\n',k);
            break

        end
        
        weight_propulsion = weight_propulsion + ...
        0.1*((thrust*0.224809 -RegConst(1))/RegConst(2)*4.44822 - weight_propulsion) ;

    end
    
    Cl_stall = airfoil_Cl_max * AR ./ (AR + 2); % finite wing correction 
    Cl_takeoff = Cl_stall/(1.1^2)+delta_Cl; % equation from 481
    wing_area_req = 2*weight_empty./(Cl_takeoff .* air_density * (Takeoff_velocity^2)); %required wing area


    err_wing_ref_area = sum(sum(abs(wing_ref_area - wing_area_req)));
    fprintf('err  %d  %f\n ', i, err_wing_ref_area);

    if err_wing_ref_area < 1e-8
        fprintf('wing_ref_area converged after %d\n',i);
        break

    end

    wing_ref_area = wing_ref_area + 0.1*(wing_area_req - wing_ref_area);  %m^2
end
%Takeoff Distance: Raymer 487
% MTOW = MTOW- (sensor * (0.141748 * g))
mu = 0.02 
K = 1./(pi*e*AR)
%Kt = thrust_to_weight-mu
Kt = thrust./MTOW-mu
Ka = air_density./2./(MTOW./wing_ref_area).*(mu.*Cl_takeoff-CD_0-K.*Cl_takeoff.^2)

takeoff_dist = (1/2/g./Ka.*log((Kt+Ka.*Takeoff_velocity.^2)./Kt))*3.28084; % answer in ft converted from m
% MTOW = MTOW+( sensor * (0.141748 * g))
%change negative to a very large number and NaN values to a very small number
MTOW(MTOW<= 0) = 1e6
MTOW(isnan(MTOW)==1) = 0.1

%change negative to a very large number and NaN values to a very small number
wing_ref_area(wing_ref_area<= 0) = inf
wing_ref_area(isnan(wing_ref_area)== 1) = 0.01


% %% ======================= Calculations ================================ %%

%Coefficients for a rough dynamic thrust curve
%T = T_2 * v^2 + T_1 * v + T_0
T_0 = thrust_to_weight*MTOW*0.7;
T_1 = -0.060;
T_2 = -0.015;

v_cruise = zeros(size(sensor));

for j = [1:size(sensor,1)]
    MTOW=MTOW-1.39
    j
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
%v_cruise(v_cruise>35) = 35;
laps_10min = (v_cruise)./(lap_length/3.28) * 60 * 10; %laps in 10 minutes
t_3laps = 3*(lap_length/3.28)./v_cruise*0.8; %time for 3 laps
%only account for well define aircraft
%for j = [1:size(pass,1)]
%    for i = [1:size(span_wing,2)]
%        if(j*3-45 > i)
%            t_3laps(j,i) = 100000000;
%            laps(j,i) = 0;
%            
%        end
%    end
%end        




M1 = zeros(size(sensor));
M1(t_3laps/60 < 5) = 1; % Mission 1 score 

M2 = (sensor./t_3laps); % mission 2 score
M2 = 1+M2./max(M2)

M3 = 2+(laps_10min .* sensorLength .* sensorWeight)./max(laps_10min .* sensorLength .* sensorWeight); % Mission 3 score
GM = 6./(sensor+5); % Ground Mission Score
score = (M1 + M2 + M3+GM);

% for i = 1:size(pass)
%     for j = 1:size(pass)
%         if span_wing(i,j)/AR(i,j)/in2m >16
%             fprintf('here')
%             score(i,j) = 0
%         end
%     end
% end

% %% ========================== Plotting ================================= %%
span_wing = span_wing/in2m
% figure(1)
% shading interp;
% [ANALY2 , ANALY2] = contourf( sensor, AR, score);
% set(ANALY2,'edgecolor','none');

% title('Scoring Analysis','FontSize',23);
% 

% ylabel('AR','FontSize',32,'FontWeight','bold');
% xlabel('# Sensors','FontSize',32,'FontWeight','bold');
% hold on 
% set(gca,'fontsize',20)

% max_score = max(max(score));
% scatter(1, 4,250, 'ko','filled')

% c = colorbar;
% c.Label.String = 'Normalized Score';
% c.Label.FontSize = 16

% fprintf('Max Score %.2f\n',max_score);
% % fprintf('RAC %.2f\n', RAC(score== max_score));
% fprintf('passengers %.2f\n', pass(score== max_score));
% fprintf('span_wing [in]%.2f\n', span_wing(score== max_score));
% fprintf('MTOW [N]%.2f\n', MTOW(score== max_score));
% fprintf('wing_ref_area [m^2]%.2f\n', wing_ref_area(score== max_score));
% fprintf('span_wingpan [in] %.2f\n', span_wing(score== max_score));
% fprintf('Chord [in] %.2f\n', span_wing(score== max_score)/AR(score== max_score));
% fprintf('AR %.2f\n', AR(score== max_score));
% fprintf('V %.2f\n', v_cruise(score== max_score));
% fprintf('Cl %.2f\n', Cl_takeoff(score== max_score));

% fprintf('RAC %.2f\n', RAC(score== max_score));
fprintf('passengers %.2f\n', sensor(1,24));
fprintf('span_wing [in]%.2f\n', span_wing(1,24));
fprintf('MTOW [N]%.2f\n', MTOW(1,24));
fprintf('EW [N]%.2f\n', MTOW(1,24)-( sensor(1,24) * (0.141748 * g)));
fprintf('wing_ref_area [m^2]%.2f\n', wing_ref_area(1,24));
fprintf('span_wingpan [in] %.2f\n', span_wing(1,24));
fprintf('Chord [in] %.2f\n', span_wing(1,24)/AR(1,24));
fprintf('AR %.2f\n', AR(1,24));
fprintf('V %.2f\n', v_cruise(1,24));
fprintf('Cl %.2f\n', Cl_takeoff(1,24)-delta_Cl);
fprintf('Thrust %.2f\n', thrust(1,24));



% figure
% contourf( sensor, span_wing, score);
% ylabel('b_s','FontSize',20);
% xlabel('# Sensors','FontSize',20);
% hold on 
% c = colorbar;
% c.Label.String = 'Normalized Score';
% c.Label.FontSize = 16

figure
contourf( sensor, span_wing, MTOW./wing_ref_area);
ylabel('b_s','FontSize',20);
xlabel('# Sensors','FontSize',20);
hold on 
c = colorbar;
c.Label.String = 'Wing Loading (MTOW/Sref)';
c.Label.FontSize = 16

figure
contourf(sensorWeight,sensorLength, score);
ylabel('Sensor Length (m)','FontSize',20);
xlabel('Sensor Weight (kg)','FontSize',20);
hold on
c = colorbar;
c.Label.String = 'Normalized Score';
c.Label.FontSize = 16

figure
contourf(M2, M3, score);
ylabel('Normalized Mission 3 Score','FontSize',20);
xlabel('Normalized Mission 2 Score','FontSize',20);
hold on 
c = colorbar;
c.Label.String = 'Total Normalized Score';
c.Label.FontSize = 16

% figure
% contourf( passengers, span_wing, (M1 + M2 + M3));
% title('(M1 + M2 + M3)','FontSize',23);
% ylabel('b_s','FontSize',20);
% xlabel('passengers','FontSize',20);
% hold on 
% c = colorbar;


% figure
% contourf( passengers, span_wing, MTOW);
% title('MTOW','FontSize',23);
% ylabel('b_s','FontSize',20);
% xlabel('passengers','FontSize',20);
% hold on 
% c = colorbar;


% figure
% contourf( passengers, span_wing, weight_empty);
% title('weight_empty','FontSize',23);
% ylabel('b_s','FontSize',20);
% xlabel('passengers','FontSize',20);
% hold on 
% c = colorbar;


% figure
% contourf( passengers, span_wing, AR);
% title('AR','FontSize',23);
% ylabel('b_s','FontSize',20);
% xlabel('passengers','FontSize',20);
% hold on 
% c = colorbar;


% figure
% contourf( passengers, span_wing, wing_ref_area);
% title('s','FontSize',23);
% ylabel('b_s','FontSize',20);
% xlabel('passengers','FontSize',20);
% hold on 
% c = colorbar;


% figure
% contourf( passengers, span_wing, v_cruise);
% title('v','FontSize',23);
% ylabel('b_s','FontSize',20);
% xlabel('passengers','FontSize',20);
% hold on 
% c = colorbar;
