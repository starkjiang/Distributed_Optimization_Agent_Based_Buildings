function [water_flow_cool_build,valve_po_heat_build,valve_po_cool_build,Actual_temp,SAT_build,MAT_build,E_actual_build,E_actual_AHU_build,E_actual_VAV_build,E_actual_fan_build,mdot_build,VP_build,delta_DAT_build,T_outside_bls_build] = Build_sim(delta_DAT,Cycle,SAT,OAT_Prog,T,mdot,SAT_sp,VP_init,MAT,valve_po_heat,valve_po_cool,T_outside_bls)
%Building simulation
%   According to the optimization result, we got the optimal supply air
%   temperature which can subsequently generate the supply air temperature
%   setpoint for the building simulation

% Update for building simulation(Every simulation would use the initial
% values of every optimization as the initial values)
global n alphac alphah c_p eta d delta_t kp_m k_I kp_dat k_I_dat kp1_dat k1_I_dat kp_dat1 k_I_dat1 xl xh Cl Ch u coef a_0 a_1 a_2 delta sigma1 sigma2 sigma3 sigma4
T((1440/Cycle)*(u-1)+1:(1440/Cycle)*(u-1)+3,:) = T(1:3,:);
SAT_AC((1440/Cycle)*(u-1)+1:(1440/Cycle)*(u-1)+3) = SAT(1:3,1);
mdot((1440/Cycle)*(u-1)+1:(1440/Cycle)*(u-1)+3,:) = mdot(1:3,:);
MAT((1440/Cycle)*(u-1)+1:(1440/Cycle)*(u-1)+3) = MAT(1:3);

T_outside_bls((1440/Cycle)*(u-1)+1:(1440/Cycle)*(u-1)+3,:) = T_outside_bls(1:3,:);
valve_position((1440/Cycle)*(u-1)+1:(1440/Cycle)*(u-1)+3,:) = VP_init(1:3,:);
delta_DAT((1440/Cycle)*(u-1)+1:(1440/Cycle)*(u-1)+3,:) = delta_DAT(1:3,:);


H_set(1:540)=65;
H_set(541:1080)=67;
H_set(1081:1443)=65;
C_set(1:540)=72;
C_set(541:1080)=70;
C_set(1081:1443)=72;



for v = (1440/Cycle)*(u-1)+3:u*(1440/Cycle)+2
    %% Dynamics between SAT and MAT (AHU)
    
    % The gross mass flow rate is the sum of every zone mass flow rate
    Mdot = mdot(v,1) + mdot(v,2) + mdot(v,3) + mdot(v,4);
%       Mdot = mdot(v,1) + mdot(v,2)+mdot(v,3)+mdot(v,4)+mdot(v,5)+mdot(v,6);
    % Mixed air temperature
%     MAT(v+1) = eta*((mdot(v,1)*T(v,1)+mdot(v,2)*T(v,2)+mdot(v,3)*T(v,3)+mdot(v,4)*T(v,4)+mdot(v,5)*T(v,5)+mdot(v,6)*T(v,6))/(mdot(v,1)+mdot(v,2)+mdot(v,3)+mdot(v,4)+mdot(v,5)+mdot(v,6)))+(1-eta)*OAT_Prog(v);
    MAT(v+1) = eta*((mdot(v,1)*T(v,1)+mdot(v,2)*T(v,2)+mdot(v,3)*T(v,3)+mdot(v,4)*T(v,4))/(mdot(v,1)+mdot(v,2)+mdot(v,3)+mdot(v,4)))+(1-eta)*OAT_Prog(v);
%       MAT(v+1) = (eta*(mdot(v,1)*T(v,1)+mdot(v,2)*T(v,2))/(mdot(v,1)+mdot(v,2)))+(1-eta)*OAT_Prog(v);
    % Cooling coil controller
    if MAT(v+1)>SAT_sp
        valve_po_heat(v+1)=0;
        valve_po_cool(v+1)=valve_po_cool(v)-kp_dat1*(SAT_AC(v)-SAT_AC(v-1))+k_I_dat1*(delta_t*((SAT_sp-SAT_AC(v))));
        
        if valve_po_cool(v+1)>1
            valve_po_cool(v+1) = 1;
        else
            if valve_po_cool(v+1)<0
                valve_po_cool(v+1)=0.01;
            end
        end
        
%         if valve_po_cool(v+1)-valve_po_cool(v)>0.03
%             valve_po_cool(v+1)=valve_po_cool(v)+0.03;
%         end
        if valve_po_cool(v+1)>0 && valve_po_cool(v+1)<xl
            water_flow_cool(v+1)=Cl/xl;
        elseif valve_po_cool(v+1)>xl && valve_po_cool(v+1)<xh
            water_flow_cool(v+1)=((Ch-Cl)/(xh-xl))*valve_po_cool(v+1)+(xh*Cl-xl*Ch)/(xh-xl);
        else
            water_flow_cool(v+1)=((1-Ch)/(1-xh))*valve_po_cool(v+1)+(Ch-xh)/(1-xh);
        end
        water_flow_heat(v+1) = 0;

        twi = 45; % cold water temperaure
        Fa = Mdot;
        Fw = water_flow_cool(v+1);
        tai = MAT(v+1);
%         options = simset('SrcWorkspace','current');
%                 sim('test',[],options)
       % sim('test');
       
% Conversion from simulink to Matlab code (1)
       SAT_AC(v+1) = SAT_AC(v)+(-30*delta*Fa-delta*sigma1*SAT_AC(v))*1/65+...
                     (-0.05*delta*Fw-delta*sigma2*SAT_AC(v))*1/65+(0.21*delta*twi-delta*sigma3*SAT_AC(v))*1/65+...
                     (0.9*delta*tai-delta*sigma4*SAT_AC(v))*1/50;
%         SAT_AC(v+1) = x.Data(end);
        
        % Valve capability
%         if  SAT_AC(v+1)-SAT_AC(v)>0.2
%             SAT_AC(v+1)=SAT_AC(v)+0.2;
%         end
    else
        % Heating coil controller
        if MAT(v+1)<SAT_sp
            valve_po_cool(v+1)=0;
            valve_po_heat(v+1)=valve_po_heat(v)-kp_dat1*(SAT_AC(v)-SAT_AC(v-1))+k_I_dat1*(delta_t*((SAT_sp-SAT_AC(v))));
            
            if valve_po_heat(v+1)>1
            valve_po_heat(v+1) = 1;
            else
              if valve_po_heat(v+1)<0
                valve_po_heat(v+1)=0.01;
            end
            end
            
%          if valve_po_heat(v+1)-valve_po_cool(v)>0.03
%             valve_po_heat(v+1)=valve_po_cool(v)+0.03;
%         end
        
            if valve_po_heat(v+1)>0 && valve_po_heat(v+1)<xl
                water_flow_heat(v+1)=Cl/xl;
            elseif valve_po_heat(v+1)>xl && valve_po_heat(v+1)<xh
                water_flow_heat(v+1)=((Ch-Cl)/(xh-xl))*valve_po_heat(v+1)+(xh*Cl-xl*Ch)/(xh-xl);
            else
                water_flow_heat(v+1)=((1-Ch)/(1-xh))*valve_po_heat(v+1)+(Ch-xh)/(1-xh);
            end
            water_flow_cool(v+1) = 0;
            
        
            twi = 120;%122 % hot water temperature
            Fa = Mdot;
            Fw = water_flow_heat(v+1);
            tai = MAT(v+1);
%             options = simset('SrcWorkspace','current');
%                 sim('test',[],options)
            %sim('test');
            
% Conversion from simulink to matlab (2)

             SAT_AC(v+1) = SAT_AC(v)+(-30*delta*Fa-delta*sigma1*SAT_AC(v))*1/65+...
                           (-0.05*delta*Fw-delta*sigma2*SAT_AC(v))*1/65+(0.21*delta*twi-delta*sigma3*SAT_AC(v))*1/65+...
                           (0.9*delta*tai-delta*sigma4*SAT_AC(v))*1/50;
%             SAT_AC(v+1) = x.Data(end);
            
        % Valve capability   
%         if  SAT_AC(v+1)-SAT_AC(v)>0.2
%             SAT_AC(v+1)=SAT_AC(v)+0.2;
%         end
            
        else
            valve_po_cool(v+1) = 0;
            valve_po_heat(v+1) = 0;
            water_flow_heat(v+1) = 0;
            water_flow_cool(v+1) = 0;
        end
    end
    
    %% Building simulation after getting the actual supply air temperature(VAV)
    for j = 1:n
        %     mdot controller  (Temperature greater than cooling setpoint)
        if T(v,j)>C_set(v);
            if mdot(v,j)>1
                mdot(v+1,j)=1;
                delta_DAT(v+1,j)=0;
                valve_position(v+1,j)=0;
            elseif mdot(v,j)<0
                mdot(v+1,j)=0;
                delta_DAT(v+1,j)=0;
                valve_position(v+1,j)=0;
            else
                mdot(v+1,j)=mdot(v,j)-kp_m(j)*(T(v,j)-T(v-1,j))+k_I(j)*(delta_t*((H_set(v)-T(v,j))));
                if mdot(v+1,j)<0
                    mdot(v+1,j)=0.01;
                else
                    if mdot(v+1,j)>1
                        mdot(v+1,j)=1;
                    end
                end
                delta_DAT(v+1,j)=0;
                valve_position(v+1,j)=0;
            end
            %    Reheat controller (Temperature less than heating setpoint)
        elseif T(v,j)<H_set(v);
            
            valve_position(v+1,j)=valve_position(v,j)-kp_dat(j)*(T(v,j)-T(v-1,j))+k_I_dat(j)*(delta_t*((H_set(v)-T(v,j))));
            
            if valve_position(v+1,j)>1
                valve_position(v+1,j)=1;
            else
                if valve_position(v+1,j)<0
                    valve_position(v+1,j)=0.01;
                end
            end
            
            if valve_position(v+1,j)>0 && valve_position(v+1,j)<xl
                water_flow(v+1,j)=Cl/xl;
            elseif valve_position(v+1,j)>xl && valve_position(v+1,j)<xh
                water_flow(v+1,j)=((Ch-Cl)/(xh-xl))*valve_position(v+1,j)+(xh*Cl-xl*Ch)/(xh-xl);
            else
                water_flow(v+1,j)=((1-Ch)/(1-xh))*valve_position(v+1,j)+(Ch-xh)/(1-xh);
            end
            twi=120;
            Fa=mdot(v,j);
            Fw=water_flow(v+1,j);
            tai=SAT_AC(v+1);
            
%             options = simset('SrcWorkspace','current');
%                 sim('test',[],options)

% Conversion from simulink to matlab
            T_outside_bls(v+1,j) = T_outside_bls(v,j)+(-30*delta*Fa-delta*sigma1*T_outside_bls(v,j))*1/65+...
                           (-0.05*delta*Fw-delta*sigma2*T_outside_bls(v,j))*1/55+(0.21*delta*twi-delta*sigma3*T_outside_bls(v,j))*1/4+...
                           (0.9*delta*tai-delta*sigma4*T_outside_bls(v,j))*1/50;
                       
%             T_outside(v+1,j)=x.Data(end);

            delta_DAT(v+1,j)=T_outside_bls(v+1,j)-SAT_AC(v+1);
            % Reheat coil capability
            if delta_DAT(v+1,j) > 20
                delta_DAT(v+1,j) = 20;
            elseif delta_DAT(v+1,j) < 0
                
                delta_DAT(v+1,j) = 0;
            end
            mdot(v+1,j)=0.25;
%         Zone temperature in between, activate the damper to cool the
%        temperature closer to heating set point
        elseif T(v,j)>H_set(v)
            if mdot(v,j)>1
                mdot(v+1,j)=1;
                delta_DAT(v+1,j)=0;
                valve_position(v+1,j)=0;
            elseif mdot(v,j)<0
                mdot(v+1,j)=0;
                delta_DAT(v+1,j)=0;
                valve_position(v+1,j)=0;
            else
                mdot(v+1,j)=mdot(v,j)-kp_m(j)*(T(v,j)-T(v-1,j))+k_I(j)*(delta_t*((H_set(v)-T(v,j))));
                if mdot(v+1,j)<0
                    mdot(v+1,j)=0.01;
                else
                    if mdot(v+1,j)>1
                        mdot(v+1,j)=1;
                    end
                end
                delta_DAT(v+1,j)=0;
                valve_position(v+1,j)=0;
            end
                
            else
                mdot(v+1,j)=0.25;
                valve_position(v+1,j)=0;
                delta_DAT(v+1,j)=0;
            end
      
        T(v+1,j) = coef(:,j)'*[T(v,j);T(v-1,j);T(v-2,j);OAT_Prog(v);mdot(v+1,j)*(SAT_AC(v+1)+delta_DAT(v+1,j)- T(v,j));d];
        %% Calculate the energy
        Eact(v+1,j) = alphac*mdot(v+1,j)*c_p*abs((SAT_AC(v+1)-MAT(v+1)))+alphah*mdot(v+1,j)*c_p*delta_DAT(v+1,j)+a_0+a_1*mdot(v+1,j)+a_2*mdot(v+1,j)^2;
        Eact_AHU(v+1,j) = alphac*mdot(v+1,j)*c_p*abs((SAT_AC(v+1)-MAT(v+1)));
        Eact_VAV(v+1,j) = alphah*mdot(v+1,j)*c_p*delta_DAT(v+1,j);
        Eact_fan(v+1,j) = a_0+a_1*mdot(v+1,j)+a_2*mdot(v+1,j)^2;
        Actual_temp(v,j)=T(v,j);
        E_actual_build (v,j) = Eact(v,j);
        E_actual_AHU_build(v,j) = Eact_AHU(v,j);
        E_actual_VAV_build(v,j) = Eact_VAV(v,j);
        E_actual_fan_build(v,j) = Eact_fan(v,j);

    end
    
    % Actual supply air temperature, supply air temperature setpoint and mixed
    % air temperature
 
   
    T_outside_bls_build(v,:) = T_outside_bls(v,:);
    SAT_build(v) = SAT_AC(v);
    MAT_build(v) = MAT(v);
    mdot_build(v,:) = mdot(v,:);
    VP_build(v,:) = valve_position(v,:);
    delta_DAT_build(v,:)=delta_DAT(v,:);
    valve_po_heat_build(v) = valve_po_heat(v);
    valve_po_cool_build(v) = valve_po_cool(v);
    water_flow_cool_build(v) = water_flow_cool(v);
end