%% Main program for the energy optimization
%%

clear all
close all
clc

%% Illustration
% SAT is the supply air temperature for optimization process whose
% dimension is 4 by 1.
% SAT_sp is the supply air temperature for building simulation which is a
% scalar.
% SAT_AC is the actual supply air temperature for AHU which is a scalar
%% Parameters

global n iter theta gamma alphac alphah alpha c_p eta d delta delta_t w S kp_m k_I kp_dat ...
                  k_I_dat kp1_dat k1_I_dat kp_dat1 k_I_dat1 xl xh Cl Ch u coef H_set C_set Pi a_0 a_1 a_2 sigma1 sigma2 sigma3 sigma4

load('Alpha_Data_ARX_1_Total.csv');
Coe1 = Alpha_Data_ARX_1_Total(:,1);
Coe2 = Alpha_Data_ARX_1_Total(:,1);
Coe3 = Alpha_Data_ARX_1_Total(:,2);
Coe4 = Alpha_Data_ARX_1_Total(:,2);

coef = [Coe1,Coe2,Coe3,Coe4];

sigma1 = 0.15;
sigma2 = 0.25;
sigma3 = 0.15;
sigma4 = 0.16;

n = 4; % the number of zones
Pi = (1/n)*ones(n); % uniform state transition matrix

Build_sim_time=15;
Cycle =1440/Build_sim_time;  % Assuming one minute sampling time
iter= 50;
theta = 0.1;
gamma = 0.1;
alphac = 0.2857; %0.010; % cooling coefficient
alphah = 1.3889; %0.020; % reheat coefficient
alpha = 1.3968; %0.024; % heating coefficient for ideal energy
c_p = 1.001; % specific heat capacity
eta = 0.85; % coefficient for mixed air temperature

%% Measurement error
d = 1;
%% Stepsize
delta_t = 1; % stepsize for controller
delta = 0.1; % stepsize for subgradient

%% Outside air temperature
load ('OAT_all.mat')
OAT_Prog = OAT_all(:,3);
%%
T_outside_opt_init = zeros(iter,n);
T_outside_opt_init(1:3,1:n) = 50;

T_outside_bls_init = zeros(1442,n);
T_outside_bls_init(1:3,1:n) = 50;


grad_SAT = zeros(iter,1);


% Supply air temperature
SAT_init(1,1:n) = 55;
SAT_init(2,1:n) = 55;
SAT_init(3,1:n) = 55;


% Weighting factor
w = 0.2;

% fan power coefficients
a_0 = 2.826;
a_1 = -2.722;
a_2 = 1.037;

% Scalar parameter
S = 100;

% Mixed air temperature
MAT_init(2:3) = 65;
% Mass flow rate
mdot_init(1:3,1:n)=0.25;

% Reheat in VAV box
delta_DAT_init(1:3,1:n)=1;

% Controller gains


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 6 zones
% kp_m=[-0.0009;0.0005;-0.0001;0.0003;-0.0004;0.0008]; k_I=[-0.0005;-0.00003;-0.0004;-0.00006;-0.0001;-0.00007];
% kp_dat=[0.5;-0.007;0.8;-0.004;0.6;-0.009]; k_I_dat=[-0.0009;0.00075;-0.0006;0.00055;-0.0002;0.00015];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4 zones
kp_m=[-0.005;-0.005;-.5;-0.005]; k_I=[-0.00009;-0.00003;-0.0009;-0.00003];
kp_dat=[-0.005;-0.0007;-0.0005;-0.0007]; k_I_dat=[-0.0009;-0.0075;-0.0009;-.00005];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%% 6 zones
% kp_dat1=0.05; k_I_dat1=0.05;
%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%% 4 zones
kp_dat1=-0.005; k_I_dat1=-0.005;
%%%%%%%%%%%%%%%

% Valve position for VAV
VP_init(1:3,1:n)=0;xl=0.3;xh=0.85;Cl=0.08;Ch=0.9907;
water_flow=[];

% xl=0.3;xh=0.85;Cl=0.08;Ch=0.9907;

% Valve position for AHU
valve_po_cool_init(1:3)=0;
water_flow_cool=[];
valve_po_heat_init(1:3)=0;
water_flow_heat=[];

% Set points for zone temperature
H_set(1:540)=65;
H_set(541:1080)=67;
H_set(1081:1443)=65;
C_set(1:540)=72;
C_set(541:1080)=70;
C_set(1081:1443)=72;


OAT_Prog(1441) = OAT_Prog(1440);
OAT_Prog(1442) = OAT_Prog(1440);
OAT_Prog(1443) = OAT_Prog(1440);


% Initial zone temperature

IC = [63;63;63;63;64;64;64;64];
T_init(1:3,1)=IC(1);T_init(1:3,2)=IC(2);
T_init(1:3,3)=IC(1);T_init(1:3,4)=IC(2);

% Outside air temperature
Zone_Temp=[];
% Building simulation to get the zone temperature, actual supply air
% temperature, mixed air temperature and actual energy
% Conducting optimization to obtain the SAT_setpoint supply air temperature plot along with iteration number
SAT_setpoint=[]; % for storage of supply air temperature set point
SAT_all_iter=[]; % for storage of supply air temperature in optimization
SAT_build_all=[];% for storage of actual supply air temperature
MAT_build_all=[];% for storage of mixed air temperature in building simulation
grad_SAT_iter_all=[]; % for storage of subgradients
E_actual_build_all = []; % for storage of actual energy consumption
E_actual_AHU_build_all = []; % for storage of actual energy consumed by AHU
E_actual_VAV_build_all = []; % for storage of actual energy consumed by VAV
E_actual_fan_build_all = []; % for storage of actual energy consumed by fan
delta_DAT_build_all = [];% for storage of reheat
valve_po_heat_build_all = []; % for storage of heating coil valve position in AHU
valve_po_cool_build_all = []; % for storage of cooling coil valve position in AHU
water_flow_cool_build_all = []; % for storage of water flow rate of cooling coil in AHU
VP_build_all = [];
mdot_build_all = [];
T_outside_bls_build_all = [];
for u = 1:Cycle          
    
    % In optimization heating and cooling set points and outside air
    % temperature are constant.
    heat_init=H_set((u-1)*15+1);
    cool_init=C_set((u-1)*15+1);
    oat_init=OAT_Prog((u-1)*15+1);
    
    % Optimization module
    [ SAT_sp,SAT_iter,grad_SAT_iter] = opt_build(heat_init,cool_init,oat_init,mdot_init,delta_DAT_init,VP_init,T_init,SAT_init,MAT_init,T_outside_opt_init);
    
    
    % Builidng simulation module
    [water_flow_cool_build,valve_po_heat_build,valve_po_cool_build,Actual_temp,SAT_build,MAT_build,E_actual_build,E_actual_AHU_build,E_actual_VAV_build,E_actual_fan_build,mdot_build,VP_build,delta_DAT_build,T_outside_bls_build] = Build_sim(delta_DAT_init,Cycle,SAT_init,OAT_Prog,T_init,mdot_init,SAT_sp,VP_init,MAT_init,valve_po_heat_init,valve_po_cool_init,T_outside_bls_init);
    % Update zone temperature, supply air temperature, mass flow rate,
    % reheat and outside air temperature for next optimization
    T_init(1:3,:) =Actual_temp(end-2:end,:);
    SAT_init(1:3,1) = SAT_build(end-2:end);SAT_init(1:3,2) = SAT_build(end-2:end);SAT_init(1:3,3) = SAT_build(end-2:end);SAT_init(1:3,4) = SAT_build(end-2:end);

    T_outside_opt_init(1:3,:) = T_outside_bls_build(end-2:end,:);
    T_outside_bls_init(1:3,:) = T_outside_bls_build(end-2:end,:);
    

    delta_DAT_init(1:3,:) = delta_DAT_build(end-2:end,:);
    mdot_init(1:3,:)=mdot_build(end-2:end,:);
    VP_init(1:3,:)=VP_build(end-2:end,:);
    MAT_init(1:3) = MAT_build(end-2:end);
    valve_po_heat_init(u*(1440/Cycle)+3) = valve_po_heat_build(end);
    valve_po_cool_init(u*(1440/Cycle)+3) = valve_po_cool_build(end);
    %%
    T_outside_bls_build_all = [T_outside_bls_build_all;T_outside_bls_build(((1440/Cycle)*(u-1)+3:u*(1440/Cycle)+2),:)];
    Zone_Temp = [Zone_Temp;Actual_temp(((1440/Cycle)*(u-1)+3:u*(1440/Cycle)+2),:)];
    SAT_build_all = [SAT_build_all;(SAT_build((1440/Cycle)*(u-1)+3:u*(1440/Cycle)+2))'];
    MAT_build_all = [MAT_build_all;(MAT_build((1440/Cycle)*(u-1)+3:u*(1440/Cycle)+2))'];
    SAT_setpoint=[SAT_setpoint;SAT_sp*ones(1440/Cycle,1)];
    SAT_all_iter=[SAT_all_iter;SAT_iter];
    grad_SAT_iter_all=[grad_SAT_iter_all;grad_SAT_iter];
    E_actual_build_all = [E_actual_build_all;E_actual_build(((1440/Cycle)*(u-1)+3:u*(1440/Cycle)+2),:)];
    E_actual_AHU_build_all = [E_actual_AHU_build_all;E_actual_AHU_build(((1440/Cycle)*(u-1)+3:u*(1440/Cycle)+2),:)];
    E_actual_VAV_build_all = [E_actual_VAV_build_all;E_actual_VAV_build(((1440/Cycle)*(u-1)+3:u*(1440/Cycle)+2),:)];
    E_actual_fan_build_all = [E_actual_fan_build_all;E_actual_fan_build(((1440/Cycle)*(u-1)+3:u*(1440/Cycle)+2),:)];
    delta_DAT_build_all = [delta_DAT_build_all;delta_DAT_build(((1440/Cycle)*(u-1)+3:u*(1440/Cycle)+2),:)];
    valve_po_heat_build_all = [valve_po_heat_build_all;(valve_po_heat_build((1440/Cycle)*(u-1)+3:u*(1440/Cycle)+2))'];
    valve_po_cool_build_all = [valve_po_cool_build_all;(valve_po_cool_build((1440/Cycle)*(u-1)+3:u*(1440/Cycle)+2))'];
    water_flow_cool_build_all = [water_flow_cool_build_all;(water_flow_cool_build((1440/Cycle)*(u-1)+3:u*(1440/Cycle)+2))'];
    VP_build_all = [VP_build_all;VP_build(((1440/Cycle)*(u-1)+3:u*(1440/Cycle)+2),:)];
    mdot_build_all = [mdot_build_all;mdot_build(((1440/Cycle)*(u-1)+3:u*(1440/Cycle)+2),:)];
end
    
% plotting
% figure
% subplot(4,1,1), plot(E_actual_build_all(:,1),'LineWidth',3);
% title('Actual energy consumed for four zones')
% legend('zone 1')
% xlabel('sampling time(min)')
% ylabel('Energy cost')
% subplot(4,1,2), plot(E_actual_build_all(:,2),'LineWidth',3);
% legend('zone 2')
% xlabel('sampling time(min)')
% ylabel('Energy cost')
% subplot(4,1,3), plot(E_actual_build_all(:,3),'LineWidth',3);
% legend('zone 3')
% xlabel('sampling time(min)')
% ylabel('Energy cost')
% subplot(4,1,4), plot(E_actual_build_all(:,4),'LineWidth',3);
% legend('zone 4')
% xlabel('sampling time(min)')
% ylabel('Energy cost')

figure
plot(Zone_Temp,'LineWidth',3)
hold on
plot(C_set,'g-','LineWidth',3)
plot(H_set,'k-','LineWidth',3)
title('Zone temperatures for four zones')
legend('zone 1','zone 2','zone 3','zone 4','zone 5','zone 6','Cooling set point','Heating set point')
xlabel('Time (h)')
ylabel('Zone temperature (F)')
ax = gca;
set(ax,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]*60);
set(ax,'XTickLabel',{'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24'});


figure
plot(SAT_all_iter,'LineWidth',3)
title('Supply air temperatures for four zones')
legend('zone 1','zone 2','zone 3','zone 4','zone 5','zone 6')
xlabel('Iteration number')
ylabel('Optimized supply air temperature')

figure
plot(MAT_build_all,'LineWidth',3)
hold on
plot(OAT_Prog,'r','LineWidth',3)
title('Mixed air temperature and Outside air temperature')
legend('Mix air Temperature','Outside air temperature')
xlabel('Time (h)')
ylabel('Temperature (F)')
ax = gca;
set(ax,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]*60);
set(ax,'XTickLabel',{'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24'});

figure
plot(SAT_build_all,'LineWidth',3)
hold on
plot(SAT_setpoint,'r','LineWidth',3)
title('Comparison of actual SAT and SAT setpoint')
legend('Actual supply air temperature','Supply air temperature setpoint')
xlabel('Time(h)')
ylabel('Supply air temperature(F)')
ax = gca;
set(ax,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]*60);
set(ax,'XTickLabel',{'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24'});

figure
plot(delta_DAT_build_all,'LineWidth',3)
legend('Zone 1','Zone 2','Zone 3','Zone 4','Zone 5','Zone 6')
title('Reheat for VAV boxes')
xlabel('Time(h)')
ylabel('Reheat for VAV boxes(F)')
ax = gca;
set(ax,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]*60);
set(ax,'XTickLabel',{'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24'});

figure
plot(mdot_build_all,'LineWidth',3)
legend('Zone 1', 'Zone 2','Zone 3', 'Zone 4','Zone 5', 'Zone 6')
title('Mass flow rate for VAV boxes')
xlabel('Time(h)')
ylabel('Mass flow rate(kg/s)')
ax = gca;
set(ax,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]*60);
set(ax,'XTickLabel',{'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24'});

% Actual energy consumption
% sum(E_actual_build_all(:,1))+sum(E_actual_build_all(:,2))+sum(E_actual_build_all(:,3))+sum(E_actual_build_all(:,4))
E_average = zeros(24,n);
E_hour_all = zeros(24,n);
for ii = 1:24
    for jj = 1:n
    E_actual = E_actual_build_all(:,jj);
    E_average(ii,jj) = mean(E_actual(1+(ii-1)*60:60*ii));
    E_hour_all(ii,jj) = 3600*E_average(ii,jj);
    end
end
% Energy_cmu=(sum(E_hour_all(:,1))+sum(E_hour_all(:,2))+sum(E_hour_all(:,3))+sum(E_hour_all(:,4))+sum(E_hour_all(:,5))+sum(E_hour_all(:,6)))/3600
Energy_cmu=(sum(E_hour_all(:,1))+sum(E_hour_all(:,2))+sum(E_hour_all(:,3))+sum(E_hour_all(:,4)))/3600

%% AHU
E_AHU_average = zeros(24,n);
E_AHU_hour_all = zeros(24,n);
for ii = 1:24
    for jj = 1:n
    E_AHU_actual = E_actual_AHU_build_all(:,jj);
    E_AHU_average(ii,jj) = mean(E_AHU_actual(1+(ii-1)*60:60*ii));
    E_AHU_hour_all(ii,jj) = 3600*E_AHU_average(ii,jj);
    end
end

% Energy_AHU=(sum(E_AHU_hour_all(:,1))+sum(E_AHU_hour_all(:,2))+sum(E_AHU_hour_all(:,3))+sum(E_AHU_hour_all(:,4))+sum(E_AHU_hour_all(:,5))+sum(E_AHU_hour_all(:,6)))/3600
Energy_AHU=(sum(E_AHU_hour_all(:,1))+sum(E_AHU_hour_all(:,2))+sum(E_AHU_hour_all(:,3))+sum(E_AHU_hour_all(:,4)))/3600

%% VAV
E_VAV_average = zeros(24,n);
E_VAV_hour_all = zeros(24,n);
for ii = 1:24
    for jj = 1:n
    E_VAV_actual = E_actual_VAV_build_all(:,jj);
    E_VAV_average(ii,jj) = mean(E_VAV_actual(1+(ii-1)*60:60*ii));
    E_VAV_hour_all(ii,jj) = 3600*E_VAV_average(ii,jj);
    end
end

% Energy_VAV=(sum(E_VAV_hour_all(:,1))+sum(E_VAV_hour_all(:,2))+sum(E_VAV_hour_all(:,3))+sum(E_VAV_hour_all(:,4))+sum(E_VAV_hour_all(:,5))+sum(E_VAV_hour_all(:,6)))/3600
Energy_VAV=(sum(E_VAV_hour_all(:,1))+sum(E_VAV_hour_all(:,2))+sum(E_VAV_hour_all(:,3))+sum(E_VAV_hour_all(:,4)))/3600

%% Fan
E_fan_average = zeros(24,n);
E_fan_hour_all = zeros(24,n);
for ii = 1:24
    for jj = 1:n
    E_fan_actual = E_actual_fan_build_all(:,jj);
    E_fan_average(ii,jj) = mean(E_fan_actual(1+(ii-1)*60:60*ii));
    E_fan_hour_all(ii,jj) = 3600*E_fan_average(ii,jj);
    end
end

% Energy_fan=(sum(E_fan_hour_all(:,1))+sum(E_fan_hour_all(:,2))+sum(E_fan_hour_all(:,3))+sum(E_fan_hour_all(:,4))+sum(E_fan_hour_all(:,5))+sum(E_fan_hour_all(:,6)))/3600
Energy_fan=(sum(E_fan_hour_all(:,1))+sum(E_fan_hour_all(:,2))+sum(E_fan_hour_all(:,3))+sum(E_fan_hour_all(:,4)))/3600

% figure
% bar(Energy_cmu,'b','stacked');hold on,bar(Energy_AHU,'r','stacked'),bar(Energy_VAV,'g','stacked'),bar(Energy_fan,'y','stacked')