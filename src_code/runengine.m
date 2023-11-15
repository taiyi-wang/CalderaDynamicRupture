function [t,x, cfail] = runengine (const, c, tspan, x0,  PeakGoal, vss, MaxComputTime)
% [t,x, cfail] = runengine (const, c, tspan, x0,  PeakGoal, vss)
%  function to run RDI solution


%  	options = odeset('RelTol',1e-7,'AbsTol',[1e-10  1e-10 1e-40 1e-10],... 
%  	'Jacobian', @calcol_Jac_RDI,'Events',@calcolEvents_down,'NonNegative',1);  
        options = odeset('RelTol',1e-6,'AbsTol',[1e-10  1e-10 1e-25 1e-10],... 
       	'Jacobian', @calcol_Jac_RDI,'Events',@calcolEvents_down);  
%    options = odeset('RelTol',1e-6,'AbsTol',[1e-10  1e-10 1e-25 1e-10]); %,... 
%     	'Jacobian', @calcol_Jac_RDI,'Events',@calcolEvents_down);  
    
    cfail = []; lastwarn(''); redtol = 0; xe = [];
    % tolerance for determining end when reaching steady state or p_out
    toler = 1e-2;
    t=0; x = x0'; smallpks = 0;
    
     odeStartTime = tic;
    % integrate until end, or v = vss, or p = pout; or lower bound.  Also
    % protect against 
    while t(end) < tspan(2) && abs(x(end,1)/const.pout -1) > toler ...
            && abs(x(end,3)/vss - 1) > toler && isempty(xe) && smallpks < 10
        % TODO: add condition that havent reached event!!
        
        sol = ode15s(@calcol_ode_RDI,[t(end) tspan(2)],x(end,:),options,c,odeStartTime,...
            MaxComputTime);
        if sol.ie == 2  % exceeded time limit
            cfail = {c; 'RunTime'};
            t = [t sol.x(2:end)]; x = [x; sol.y(:,2:end)']; 
            return
        end
        
        xe = sol.xe;
        [~, msgid] = lastwarn;
        if strcmp(msgid,'MATLAB:ode15s:IntegrationTolNotMet')
            %disp('switching to ode23s')
            %cfail = 23;
            lastwarn('');
            %options = odeset(options,'RelTol',1e-4);
            %odeStartTime = tic;
            sol = ode23s(@calcol_ode_RDI,[t(end) tspan(2)],x(end,:),options,c,odeStartTime,...
            MaxComputTime);
            if sol.ie == 2  % exceeded time limit
                cfail = {c; 'RunTime'};
                t = [t sol.x(2:end)]; x = [x; sol.y(:,2:end)']; 
                return
            end
            xe = sol.xe;
            [~, msgid] = lastwarn;
            if strcmp(msgid,'MATLAB:nearlySingularMatrix')
                cfail = {c; 'nearsingular'}; 
            end
            if strcmp(msgid,'MATLAB:ode23s:IntegrationTolNotMet')
                t = [t sol.x(2:end)];
                x = [x; sol.y(:,2:end)'];  
                cfail = {c;'23intnotmet'};
                return
            end          
        end

    t = [t sol.x(2:end)];
    x = [x; sol.y(:,2:end)'];  
    %figure; semilogy(t,x(:,3)); drawnow; hold on
    [pks, pkloc, ~] = findpeaks(x(:,3),t,'MinPeakHeight',0.1);
    % catch peaks that dont become inertial
    smallpks = length(findpeaks(x(:,3)))-length(pks);
        if length(pks) > 1
            tspan(2) = tspan(2) + (PeakGoal - length(pks))*mean(diff(pkloc)) ;
        elseif length(pks) < 2 && redtol ==0 
            tspan(2) = 1.5*PeakGoal*tspan(2);
        end
    end

end
function out = calcol_ode_RDI (t,x,c,e1,e2)
%       out = calcol_ode_I (t,x,c)
%
%       t = times
%       x = model vector
%	x(1) is pressure
%	x(2) is state  variable
%	x(3) velocity
%   x(4) is slip
%
%   Version with both intertia and radiation damping
%   Recall these are all non-dimensional
% TODO: put in a flag for aging vs slip law

   if ~isreal(x)
       keyboard
   end

   a = c(1); b = c(2); d_c = c(4); sigma = c(5); alpha = c(6); 
  f0 = c(3);  kappa = c(8); gamma = c(9); pout = c(10); 
  V0 = c(12);  Gamma = c(11);  eta = c(7); astar = c(13);

    %keyboard
    p     = x(1);
	theta = x(2);
    v     = x(3);
    u     = x(4);

      
     %return rates: pressure
     out(1) =  -Gamma*(p - pout - gamma*u) + kappa*v;

     % state evolution
     if (c(14)  == 1)  % aging  law
         out(2) =  1 - theta*v/d_c;
     elseif  (c(14)  == 2)  % slip law    
         out(2) =   -theta*v/d_c*log(theta*v/d_c);
     end
 
     % ode form, unregularized    
     %out(3) =  -1./(eta + sigma*a/v).*(alpha*out(1) + sigma*b*out(2)/theta);
     %tau = sigma*(f0 + a*log(v) + b*log(theta/d_c) );
     
     % stress, asinh form
     %tau1 = a*sigma*asinh(0.5*v*exp(f0/a + b*log(theta/d_c)/a ));
     tau = a*sigma*asinh(0.5*v*exp(f0/a)*(theta/d_c)^(b/a));

     % velocity
     out(3) =  (alpha - alpha*p - tau - eta*v)/astar;
     out(4) = v;

     out = out';
end


function J = calcol_Jac_RDI (t,x,c,e1,e2)
%       out = calcol_Jac_RDI (t,x,c)
%
%       t = times
%       x = model vector
%	x(1) is pressure
%	x(2) is state  variable
%	x(3) velocity
%   x(4) is slip
%
%   Version with  intertia and radiation damping
%   Recall these are all non-dimensional
% TODO: put in a flag for aging vs slip law



   a = c(1); b = c(2); d_c = c(4); sigma = c(5); alpha = c(6); 
  f0 = c(3);  kappa = c(8); gamma = c(9);  
  Gamma = c(11);   eta = c(7); %pout = c(10);V0 = c(12);  
  astar = c(13);

    
    p     = x(1);
	theta = x(2);
    v     = x(3);
    u     = x(4);
    
    J = zeros(4,4);
    J(1,:) = [-Gamma, 0, kappa, Gamma*gamma];
    % state evolution
     if (c(14)  == 1)  % aging  law    
         J(2,:) = [0, -v/d_c, -theta/d_c, 0];
     elseif (c(14)  == 2)  % slip law   
         J(2,:) = [0, -v/d_c*(log(theta*v/d_c) + 1), -theta/d_c*(log(theta*v/d_c) + 1), 0];
     end
    
	%J(3,:) = [-alpha/astar, -sigma*b/(theta*astar), -sigma*a/(v*astar), 0];
    
    % for asinh form
    %z = 0.5*v*exp(f0/a + b*log(theta/d_c)/a);
    z = 0.5*v*exp(f0/a)*(theta/d_c)^(b/a);

    zterm = (z^2 + 1)^(-0.5);
    J(3,:) = [-alpha/astar, -a*sigma/astar*zterm*z*b/a/theta, ...
                            (-a*sigma/astar*zterm*z/v -eta/astar),  0];
    
 	J(4,:) = [0, 0, 1, 0];
  
end