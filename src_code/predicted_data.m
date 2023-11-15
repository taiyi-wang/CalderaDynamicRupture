function [dhat, cfail, varargout] = predicted_data (m, const, doplot,force,etafact)
% function dhat = predicted_data (m, const, tspan, doplot,force)
% subroutine to compute predicted data given a model;
% Modify Oct 5, 2020: make b-a the variable, call new predicted data
% Modified Dec, 2020 to remove radiation damping; const.eta = 0;
% Modified Jan 1, 2021 to make most priors linear, except d_c
% Modified Jan 6, 2021 to make stress "observations" in  MPa

radFlag = 1; % default is with radiation damping
if nargin == 2 
   doplot = 0;
end
if nargin ==3
   doplot = 0;
   radFlag = 0; % 0 for no radiation damping
end 
if nargin < 4
    force = 0;
    etafact = 1;
end

    %convert models to const structure for input into prediction 
    const = models2const(m, const);
    
    % the following lines were added specifically for deomonstrating
    % inversion with/without radiation damping
    if radFlag == 0
        const.eta = 0; % no radiation damping
    end
    if radFlag == 1
        %const.pstar = const.pstar*2;
    end 
    

    const.eta = etafact*const.eta;

    % output large value if no solution or steady-state does not exist
    Large = 1e10;
    
%  test if steady-state exists
    epsilon = 1e-5;
    foo = logspace(-40,10,1000);
    C = const.alpha*(1-const.pout) - const.sigma*const.f0;
    D = const.sigma*(const.B-const.A);  
    y = C-const.alpha*const.kappa*foo/const.Gamma + D*log(foo);
    if max(y) <  -epsilon && force ~= 1
       dhat = [Large; Large; Large; Large; Large; Large; Large];
       %cfail = {c; 'No SS'};
       cfail = [];
      return
    end

%  test if system is linearly unstable
    [~, vss, ~, ~, ~] = steadystate(const,1);
    if (const.kappa > kappac (const, vss) || const.kappa < 0 ) && force ~= 1
       dhat = 2*[Large; Large; Large; Large; Large; Large; Large];
       %cfail = {c; 'Linearlystable'};
       cfail = [];
       return
    end
    
    [c, tspan, x0] = setup (const);    
% check that initial pressure exceeds pout
    if (x0(1) < const.pout || const.pout < 0) && force ~= 1
       dhat = 3*[Large; Large; Large; Large; Large; Large; Large];
       %cfail = {c; 'p0TooSmall'};
       cfail = [];
       return        
    end
    
% Only if passes these tests compute solution
    PeakGoal = 5;  % get at least 4 peaks in v
    MaxComputTime = 30;  % maximum integration time in seconds

    
    [t,x, cfail] = runengine(const, c, tspan, x0, PeakGoal, vss, MaxComputTime);
    
    
    CharDim = [const.L,const.tstar*86400, const.pstar];  % characteristic variables in MKS units
    [T,Te,DeltaP,slip_med,Tp,tau_peak,tau_min] = prediction_dim(t, x, const, CharDim);
    
    if (isnan(T) && force ~= 1)
       %cfail = {c; 'NaN Period'};
       cfail = [];
       dhat = 5*[Large; Large; Large; Large; Large; Large; Large];
       return
    end

    % Predicted data
    % units [days, seconds, meters, MPa, days]
    dhat = [T; Te;  DeltaP; slip_med; Tp; tau_peak; tau_min];  
   
    % for plotting dimensional time series output
    if(doplot)
    disp(num2str(dhat))
    calcol('plotout_dim', const,t,x,const.pstar,const.L,const.tstar)
    end
    
    if  nargout == 4
        varargout{1} = t;
        varargout{2} = x;
    end
    
    if  nargout == 5
        varargout{1} = t;
        varargout{2} = x;
        varargout{3} = const;
    end
        
end


function [T,Te,DeltaP,slip_med,Tp,tau_peak,tau_min] = prediction_dim (t, x, const, CharDim)
%[T,Te,DeltaP,slip_med,Tp,tau_peak,tau_min] = prediction_dim (t, x, const, CharDim)
%
% function to pull out predictions, dimensionalize FIRST
% Input:
%       t, x, solution from ODE solver
%       const structure with parameters
%       CharDim vector of characteristic dimensions (L*, t*, p*)
% Output:
%       T:      median repeat time  [days]
%       Te:     measure of event duration [seconds]
%       slip:   median slip in events [meters]
%       DeltaP  median pressure increase [Mpa]
%       Tp      estimate of pressure decay [days]
%
%  P Segall November 2019
    if length(t) < 100
        T =1e10; Te = 1e10; DeltaP = 1e10; slip_med = 1e10; Tp = 1e10;
        tau_peak=1e10; tau_min=1e10;
        return
    end

    % first redimensionalize
    v = x(:,3)*CharDim(1)/CharDim(2);
    p = x(:,1)*CharDim(3);

    % to find peaks in velocity; minimum peak height is 0.1 m/s, 
    [~, pkloc, ~] = findpeaks(v,t,'MinPeakHeight',0.1);
    
    % if not more than 3 events return NaN
    if length(pkloc) < 3 %disp('warning: 3 peaks not detected');
        T =1e10; Te = 1e10; DeltaP = 1e10; slip_med = 1e10; Tp = 1e10;
        tau_peak=1e10; tau_min=1e10;
        return
    end
    
    % find median repeat time
    T = CharDim(2)/86400*median(  diff(pkloc(2:end)) );
    
    
    %To get event duration
    I = find(v >= 0.1); % all indices with  velocities > 0.1 m/s
    % jumps occur when I increments by more than 1
    II = find(diff(I)>1);
    % duration of high slip episodes is thus
    dt(1) = t(I(II(1)))-t(I(1));
    for k = 2:length(II)
        dt(k) = t(I(II(k))) - t(I(II(k-1)+1));
    end
    Te = CharDim(2)*median(dt);  % to convert to dimensional in seconds
    
    
%         % for plotting
%         figure;
%         semilogy(I,v(I),'.'); hold on
%         semilogy(I(II),v(I(II)), 'o')  % these points mark end of high velocity
%         figure;
%         semilogy(t(I),v(I),'.'); hold on
%         semilogy(t(I(II)),v(I(II)), 'o')  % these points mark end of high velocity

    
    % To get slip and  pressure change...First to get indices of peaks
        [~, pkI] = findpeaks(v,'MinPeakHeight',0.1);
        slip_med = CharDim(1)*median(diff(x(pkI,4)));
        
    % pressure;  find peaks in pressure
        [pkp] = median(findpeaks(p));  %maxima in  p
        [pkm] = median(-findpeaks(-p)); %minimum in  p
        DeltaP = 1e-6*(pkp - pkm);
 
%         [foo, pploc1, ~] = findpeaks(p,t)
%         [foo, pploc2, ~] = findpeaks(-p,t)
%         % new estimate of duration is
%         Te2 = (pploc1(end)-pploc2(end))*CharDim(2);
       
        
    % stress drop
        %tau = const.sigma*(const.f0 + const.A*log(x(:,3)) + ...
        %const.B*log(x(:,2)/const.d_c));
    
        %tau = a*sigma*asinh(0.5*v*exp(f0/a)*(theta/d_c)^(b/a));
        % NOTE: terms inside the asinh are non-dimensional, which is ok
        % since the normalization constant v0 is absorbed there
        tau = const.A*const.sigma*asinh(0.5*x(:,3).*exp(const.f0/const.A)...
            .*(x(:,2)/const.d_c).^(const.B/const.A));
        
        tau_range = max(tau) - min(tau);
        tau_peak = 1e-6*CharDim(3)*median(findpeaks(tau,'MinPeakProminence',tau_range/2));
        tau_min = 1e-6*CharDim(3)*median(-findpeaks(-tau,'MinPeakProminence',tau_range/2));
%        Delta_tau = (tau_peak - tau_min);
 
	%To get the characteristic decay time in pressure
    % indices of last cycle
    
    pp =  p(pkI(end-1):pkI(end)); tp = t(pkI(end-1):pkI(end));
    [~, Ipmax] = max(pp);  % find max and keep only decay part
    pp = pp(Ipmax:end); tp = tp(Ipmax:end);
    % redimensionalize time to days
    tp = CharDim(2)*(tp - tp(1))/86400;
    % robust fit to exponential
    %frob = robustfit(tp, log(pp))
%     clear II
     pp = pp - min(pp); pp = pp/pp(1);
%     II = find(tp < 0.8*T);
%     bls = regress(log(pp(II)),[ones(length(tp(II)),1) tp(II)']);
%     Tp = -1/bls(2);
    fun = @(x,xdata)x(1)*exp(-xdata/x(2)) + x(3);
    x0 = [1,0.5,0];
    opts = optimoptions(@lsqcurvefit,'FunctionTolerance',1e-5,'Display','off');
    [x,~,~,EXITFLAG] = lsqcurvefit(fun,x0,tp',pp,[],[],opts);
    % should be in something for exitflag, ~= 1 or 3
    Tp = x(2);
    %EXITFLAG  
    % dimensionalize Tp to days
    %Tp = CharDim(2)/86400*Tp; already in units of days 
    %keyboard
    return
    figure; plot(tp, pp, 'LineWidth',3); hold on
    plot(tp',fun(x,tp'),'r--','LineWidth',2)
    plot(tp, exp(bls(1))*exp(bls(2)*tp))
end

function signom = nominal_sig(rho_c, L)
% nominal average normal stress
%
% constants:
    rho_w = 1e3;
        g = 9.8;
      h_w = 500;  %water table depth

signom = 0.5*(rho_c-rho_w)*g*L + 1e3*g*h_w;
end

function const = models2const (m, const)
% const = models2const (m, const);
% Input models m for MCMC
% Output  existing const structure (dimensionless variabels)

% modified 10-5-20 to make b-a parameter
% modified 12-24-20 to set eta = 0; no radiation damping%
% Modified 1-1-21 to to to mostly linear parameters

% dimensionless parameters
    const.A = m(1);  const.f0 = m(4);
    const.B = m(2) + const.A;

%  characteristic variables 
    tstar = m(9); % time scale in days
    L  = m(6); %1e3; %length scale in m 
    
% for these parameters we are going to use nominal values; all units MKS
    rho = 2.5e3; g = 9.8; rho_c = 2.5e3; R = 0.8*1e3;
    pstar = rho_c*g*L;   % pressure in Pa
    signom = nominal_sig(rho_c, L);
    vs = 3.3e3;         % shear wave speed in m/s
    const.eta = 0.5*vs/g/(tstar*86400);
 
    
% Other Non-dimensional parameters to pass to solver
    const.alpha = R/2/L;
    const.d_c = 10^(m(3))/L;    % d_c in meters -> Non-dim
    const.sigma = m(5)*signom/pstar;   % sigma in Pa -> Non-dim m(5) is multiplicative factor
    const.signom = signom;
    %const.V0 = 10^(m(8))/L^3;   % V0 in cubic meters  -> Non-dim
    const.V0 = m(8)/L^3;   % V0 in cubic meters  -> Non-dim
    
    % beta m(7) is in units of l( Pa^-1)
    const.kappa = pi*R^2/(rho*g*m(8)*m(7)); 
    
    % Note pout determined uniquely by L and the 1 km height differenc 
    % between summit and eruption site
    
   delh = 0.8*1e3;
   const.pout = 1 -delh/L;
    
    % add characteristic dimesions to const structure
    const.L = L; const.tstar = tstar; const.pstar = pstar;
    
    % For full inertia include astar
    const.astar = 0.5*R/g/(tstar*86400)^2;

    
end



function c = const2vec (const)
%convert const stucture to vector
c(1) = const.A; c(2) = const.B; c(3) = const.f0; c(4) = const.d_c;
c(5) = const.sigma; c(6) = const.alpha; c(7) = const.eta; c(8) = const.kappa;
c(9) = const.gamma; c(10) = const.pout; c(11) = const.Gamma;  c(12) = const.V0;
c(13) = const.astar;
c(14) = const.evol;
%const.evol = 1;   % state evolution law: 1 = aging, 2 = slip
end

function const = cvec2const(c, const)
const.A = c(1); const.B = c(2); const.f0 = c(3);  const.d_c = c(4);
const.sigma = c(5); const.alpha = c(6); const.eta = c(7); const.kappa = c(8);
const.gamma = c(9);  const.pout = c(10); const.Gamma = c(11);  const.V0 = c(12);
const.astar = c(13);
const.evol = c(14);
end

function [Nroots, vss, pss, vss2, pss2] = steadystate (const,noout)
    % Subroutine to find steady-state solution. There can be either 0,1,or
    % 2 roots. That is passed back as Nroots.  (vss,pss) and (vss2, pss2)
    % are the steady state velocity and pressure if there are 1 or 2 roots
    % nout = true if supress printing
    
    if nargin == 1; noout = false; end
    
    epsilon = 1e-5;  % threshold for determining number of roots; could be improved
    
    % define constants
    d_c = const.d_c; A = const.A; B = const.B; sigma = const.sigma; f0 = const.f0;
    kappa = const.kappa; alpha = const.alpha; eta = const.eta; pout = const.pout;   
    
    C = alpha*(1-pout)-sigma*f0;   
    options = optimset('TolX',1e-12,'TolFun',1e-8);
    [vmax,minval,EXITFLAG] = fminbnd(@(x) -(C-(alpha*kappa +eta)*x+sigma*(B-A)*log(x)),1e-12,1e-2,...
            options);
    if EXITFLAG ~= 1
        disp('exitflag in fminbnd = '); EXITFLAG
    end
        
        
    if ~noout
    %[vmax,minval,EXITFLAG] = fminbnd(@(x) -(C-alpha*kappa*x+sigma*(B-A)*log(x)),1e-12,1e2);
    disp(['Max Value:  ',num2str(-minval)])
    disp('If > 0 2 roots, < 0 no roots, = 0 one root')
    %semilogx(vmax, -minval, 'r*')
    end
    if -minval < -epsilon
        Nroots = 0; vss = NaN; pss = NaN; vss2 = NaN; pss2 = NaN;
        return
    end

    % Handle case with multiple roots
    if -minval > epsilon
        Nroots = 2;
        
        %take larger root; these bounds seem to work
        bounds = [vmax, 1e5];
        % to experiment with smaller root
        %bounds = [1e-12, vmax];
        
        [vss, pss] = sssol(const, bounds);
        [vss2, pss2] = sssol(const, [1e-30, vmax]);


    elseif abs(minval) < epsilon % single root
        Nroots = 1;
        disp('single root')
        vss = vmax;
        pss = pout + kappa*vss;
        vss2 = []; pss2 = [];
    end
    
    if Nroots > 0 && ~noout
        d_c_crit = sigma*(B-A)/(  alpha*kappa + sigma*A/vss);
        % smaller d_c causes instability; larger causes stability
         if d_c < d_c_crit
            disp('linearly unstable')
         else
             disp('linearly stable')
         end
    end
%         vss = []; pss = [];
%     end
         
end

function [vss, pss] = sssol (const, bounds)
    % find the steady-state solution given some bounds in velocity
    % define constants
    A = const.A; B = const.B; sigma = const.sigma; f0 = const.f0;
    kappa = const.kappa; alpha = const.alpha;  pout = const.pout;  
    
  
    C = alpha*(1-pout) - sigma*f0;
    M = alpha*kappa;
    D = sigma*(B-A);
    
    if sign(C-M*bounds(1)+D*log(bounds(1))) == sign(C-M*bounds(2)+D*log(bounds(2)))
        vss = NaN; pss = NaN;
        return
    else
  
    [vss,~, ~] = fzero(@(x) C-M*x+D*log(x),bounds);
    end
    pss = kappa*vss + pout;
    
        if(0)  % if fzero fails because no root is found plot to see what happend
            vs = logspace(log10(bounds(1)), 1, 100);
            figure
            semilogx(vs,  alpha*(1-pout) - alpha*kappa*vs, 'LineWidth', 2)
             hold on
            semilogx(vs, sigma*f0 + sigma*(A-B)*log(vs), 'LineWidth', 2)
            semilogx(vs, sigma*(B-A)*log(vs) + alpha*(1-pout)-sigma*f0 - alpha*kappa*vs, ...
            'k','LineWidth', 2)
            legend('Load - Pressure','Friction', 'Total')
            ylim([-1 1])
        end
end

function [kappa_c, omega_c] = kappac (const, vss)
% kappa_c = kappac (const, vss)
% compute critical value of kappa

if vss <= 0
    disp('No steady state solution')
    return
end

kappa_c = const.sigma/const.alpha*((const.B-const.A)/const.d_c - ...
    const.A*const.Gamma/vss);

omega_c = sqrt( (const.B-const.A)/const.A - const.B/const.A*const.Gamma* ...
    const.d_c/vss)*vss/const.d_c;

end


function out = calcol_ode (t,x,c)
%       out = calcol_ode(t,x,c)
%
%       t = times
%       x = model vector
%	x(1) is pressure
%	x(2) is state  variable
%	x(3) velocity
%   x(4) is slip

% TODO: use flag to determine which state evolution law to use
%

%  unwrap constantsl this version uses vector c not structure
%c(1) = const.A; c(2) = const.B; c(3) = const.f0; c(4) = const.d_c;
%c(5) = const.sigma; c(6) = const.alpha; c(7) = const.eta; c(8) = const.kappa;
%c(9) = const.gamma; c(10) = const.pout; c(11) = const.Gamma;  c(12) = const.V0;

   a = c(1); b = c(2); d_c = c(4); sigma = c(5); alpha = c(6); 
   eta = c(7); kappa = c(8); gamma = c(9); pout = c(10); 
   Gamma = c(11);%V0 = c(12); f0 = c(3); % 
    
    p     = x(1);
	theta = x(2);
    v     = x(3);
    u     = x(4);
    

     %return rates
     %k = calck(u,uc,sig,knom);
    % k = 1;
     
     % TESTING ONLY
%      dV = 4*pi*alpha^2*u;
%      k = V0/(V0-dV);
     out(1) =  -Gamma*(p - pout - gamma*u) + kappa*v;
     
     % new form that accounts for changing chamber volume
     %out(1) =  out(1)./(1 - 4*pi*alpha^2*u/V0);

     % state evolution
      out(2) =  1 - theta*v/d_c;
     %slip law
     % out(2) =   -theta*v/d_c*log(theta*v/d_c);
     
     
     % TEST: experiment with changing normal stress
     % FAIL: this seems to require an DAE solver to really perserve
     % momentum balance. set constant to 0.
     %sigma = const.sigma + const.sigUfact*u;
     
     % ode form, unregularized    
     out(3) =  -1./(eta + sigma*a/v).*(alpha*out(1) + sigma*b*out(2)/theta);
     out(4) = v;

     out = out';
     return
end




function [c, tspan, x0] = setup (const)
    % set parameters: first as structure then convert to vector. Legacy
    c = const2vec(const);
    tspan = [0 1];
  
%     if(0)   % this was the way I used to initialize    
%         v0      = 1e-8;
%         p0      = 0.6;
%         tau0 = const.alpha*(1-p0) - const.eta*v0;    %  dv/dt = 0
%         theta0  = const.d_c * exp (  (tau0/const.sigma - const.f0 - const.A*log(v0) )/const.B );
%         %foo = sinh(tau0/(const.A*const.sigma))/(0.5*v0*exp(const.f0/const.A));
%         %theta0 = const.d_c*foo^(const.A/const.B);
%     else
        % choose low velocity and make sure tau is above steady-state.
        v0 = 1e-4;
        tau0 = 1.001*const.sigma*(const.f0 + (const.A-const.B)*log(v0) );
        theta0  = const.d_c * exp (  (tau0/const.sigma - const.f0 - const.A*log(v0) )/const.B );
        % better check that p0 is above pout
        p0 = 1 - tau0/const.alpha - const.eta*v0/const.alpha;
        mb = const.alpha*(1-p0) - tau0- const.eta*v0;
        if abs(mb) > 5e-6; disp('no momentum balance'); keyboard; end
        %foo =  p0-const.pout;
        %if foo < 0; disp('initial pressure < outlet'); end
        
%   end    
    x0      =   [p0; theta0; v0; 0];

end

