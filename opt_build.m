function [ SAT_sp,SAT_iter,grad_SAT_iter] = opt_build(h_set,c_set,oat,mdot,delta_DAT,valve_position,T,SAT_iter,MAT,T_outside_opt)
%Optimizing supply air temperature
%   Using generalzied gossip algorithm to calculate the optimal supply 
%   temperature
global  n iter theta gamma alphac alphah alpha c_p eta d delta delta_t w S kp_m k_I kp_dat k_I_dat kp1_dat k1_I_dat xl xh Cl Ch coef H_set C_set Pi a_0 a_1 a_2 sigma1 sigma2 sigma3 sigma4

T_outside_opt(1:3,:) = T_outside_opt(1:3,:);

for i = 3:iter-1
    
    
    for j = 1:n
        
        %     mdot controller  (Temperature greater than cooling setpoint)
        if T(i,j)>c_set;
            if mdot(i,j)>1
                mdot(i+1,j)=1;
                delta_DAT(i+1,j)=0;
                valve_position(i+1,j)=0;
            elseif mdot(i,j)<0
                mdot(i+1,j)=0;
                delta_DAT(i+1,j)=0;
                valve_position(i+1,j)=0;
            else
                mdot(i+1,j)=mdot(i,j)-kp_m(j)*(T(i,j)-T(i-1,j))+k_I(j)*(delta_t*((h_set-T(i,j))));
                if mdot(i+1,j)<0
                    mdot(i+1,j)=0;
                else
                    if mdot(i+1,j)>1
                        mdot(i+1,j)=1;
                    end
                end
                delta_DAT(i+1,j)=0;
                valve_position(i+1,j)=0;
            end
            %    Reheat controller (Temperature less than heating setpoint)
        elseif T(i,j)<h_set;
            
            valve_position(i+1,j)=valve_position(i,j)-kp_dat(j)*(T(i,j)-T(i-1,j))+k_I_dat(j)*(delta_t*((h_set-T(i,j))));
            
            if valve_position(i+1,j)>1
                valve_position(i+1,j)=1;
            else
                if valve_position(i+1,j)<0
                    valve_position(i+1,j)=0.01;
                end
            end
            
            if valve_position(i+1,j)>0 && valve_position(i+1,j)<xl
                water_flow(i+1,j)=Cl/xl;
            elseif valve_position(i+1,j)>xl && valve_position(i+1,j)<xh
                water_flow(i+1,j)=((Ch-Cl)/(xh-xl))*valve_position(i+1,j)+(xh*Cl-xl*Ch)/(xh-xl);
            else
                water_flow(i+1,j)=((1-Ch)/(1-xh))*valve_position(i+1,j)+(Ch-xh)/(1-xh);
            end
            
            twi=120;
            Fa=mdot(i,j);
            Fw=water_flow(i+1,j);
            tai=SAT_iter(i,j);
%             options = simset('SrcWorkspace','current');
%                 sim('test',[],options)
            %sim('test');
            
% Conversion from simulink to matlab (Assign a weight to each input and
% then convert them to discrete-time domain model. Feel free to change if
% necessary
            T_outside_opt(i+1,j) = T_outside_opt(i,j)+(-30*delta*Fa-delta*sigma1*T_outside_opt(i,j))*1/65+...
                           (-0.05*delta*Fw-delta*sigma2*T_outside_opt(i,j))*1/65+(0.21*delta*twi-delta*sigma3*T_outside_opt(i,j))*1/65+...
                           (0.9*delta*tai-delta*sigma4*T_outside_opt(i,j))*1/50;

%             T_outside_opt(i+1,j)=x.Data(end);
           
            delta_DAT(i+1,j)=T_outside_opt(i+1,j)-SAT_iter(i,j);
            % reheat coil capability
            if delta_DAT(i+1,j) > 20
                delta_DAT(i+1,j) = 20;
            else 
                delta_DAT(i+1,j) < 0
                delta_DAT(i+1,j) = 0;
            end
            mdot(i+1,j)=0.25;
            
        else
            % Make zone temperature closer to heating set point
            if T(i,j)>h_set
             
            if mdot(i,j)>1
                mdot(i+1,j)=1;
                delta_DAT(i+1,j)=0;
                valve_position(i+1,j)=0;
            elseif mdot(i,j)<0
                mdot(i+1,j)=0;
                delta_DAT(i+1,j)=0;
                valve_position(i+1,j)=0;
            else
                mdot(i+1,j)=mdot(i,j)-kp_m(j)*(T(i,j)-T(i-1,j))+k_I(j)*(delta_t*((h_set-T(i,j))));
                if mdot(i+1,j)<0
                    mdot(i+1,j)=0;
                else
                    if mdot(i+1,j)>1
                        mdot(i+1,j)=1;
                    end
                end
                
                delta_DAT(i+1,j)=0;
                valve_position(i+1,j)=0;
            end
           
                
            else
                mdot(i+1,j)=0.25;
                valve_position(i+1,j)=0;
                delta_DAT(i+1,j)=0;
            end
        end
        % Zone temperature
        T(i+1,j) = coef(:,j)'*[T(i,j);T(i-1,j);T(i-2,j);oat;mdot(i+1,j)*(SAT_iter(i,j)+delta_DAT(i+1,j)- T(i,j));d];
        % Mixed air temperature
%         MAT(i+1) = (eta*(mdot(i,1)*T(i,1)+mdot(i,2)*T(i,2)+mdot(i,3)*T(i,3)+mdot(i,4)*T(i,4)+mdot(i,5)*T(i,5)+mdot(i,6)*T(i,6))/(mdot(i,1)+mdot(i,2)+mdot(i,3)+mdot(i,4)+mdot(i,5)+mdot(i,6)))+(1-eta)*oat;
          MAT(i+1) = (eta*(mdot(i,1)*T(i,1)+mdot(i,2)*T(i,2)+mdot(i,3)*T(i,3)+mdot(i,4)*T(i,4))/(mdot(i,1)+mdot(i,2)+mdot(i,3)+mdot(i,4)))+(1-eta)*oat;
%           MAT(i+1) = (eta*(mdot(i,1)*T(i,1)+mdot(i,2)*T(i,2)))/(mdot(i,1)+mdot(i,2))+(1-eta)*oat;
        % Zone comfort cost
%         C(i+1,j) = w*(T(i+1,j)-(C_set(i+1)+H_set(i+1))/2)^2;
%         C(i-1,j) = w*(T(i-1,j)-(C_set(i-1)+H_set(i-1))/2)^2;
        
        C(i+1,j) = S*w*(T(i+1,j)-H_set(i+1))^2;
        C(i-1,j) = S*w*(T(i-1,j)-H_set(i-1))^2;
        % Waste energy cost
%         E(i+1,j) = (1-w)*(mdot(i+1,j)*c_p*((alpha-alphac)*MAT(i+1)+(alphac-alpha)*SAT_iter(i,j)+(alphah-alpha)*delta_DAT(i+1,j)))^2;
%         E(i-1,j) = (1-w)*(mdot(i-1,j)*c_p*((alpha-alphac)*MAT(i-1)+(alphac-alpha)*SAT_iter(i-2,j)+(alphah-alpha)*delta_DAT(i-1,j)))^2;
        % Total energy cost
        E(i+1,j) = alphac*mdot(i+1,j)*c_p*abs((SAT_iter(i,j)-MAT(i+1)))+alphah*mdot(i+1,j)*c_p*delta_DAT(i+1,j)+a_0+a_1*mdot(i+1,j)+a_2*mdot(i+1,j)^2;
        E(i-1,j) = alphac*mdot(i-1,j)*c_p*abs((SAT_iter(i-2,j)-MAT(i-1)))+alphah*mdot(i-1,j)*c_p*delta_DAT(i-1,j)+a_0+a_1*mdot(i-1,j)+a_2*mdot(i-1,j)^2;
        % sub-gradient
        % Central difference
        grad_SAT(i,j) = ((E(i+1,j)+C(i+1,j))-(E(i-1,j)+C(i-1,j)))/(2*delta);
        
        % Four zone information of supply air temperature and subgradients
  
        grad_SAT_iter(i,j) = grad_SAT(i,j);
        
    end
    
    % Supply air temperature update
    
    
    SAT_int_iter(:,i+1) = (1-theta)*Pi*SAT_iter(i,:)' + theta*(SAT_iter(i,:)' - gamma*grad_SAT_iter(i,:)');
    SAT_iter(i+1,:)=SAT_int_iter(:,i+1)';
   
    
    % Upper and lower bounds for supply air temperature
    for q = 1:n
        if SAT_iter(i+1,q)>80
            SAT_iter(i+1,q) = 80;
        elseif SAT_iter(i+1,q)<45
            SAT_iter(i+1,q) = 45;
        end
    end   
end

SAT_sp = mean(SAT_iter(end,:));