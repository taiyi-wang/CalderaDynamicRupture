%% Demonstrate the effect of radiation damping on collapse dynamics via a modified lumped model
% The provided time domain solution for slip is obtained by solving Eqn.C5

current_dir  = pwd;
idcs   = strfind(current_dir ,'/');
above_dir = current_dir(1:idcs(end)-1); % one level above current directory

%% set up model parameters
asp_0 = 0.5;          % default value of chamber aspect ratio
epsilon_0 = 1;        % default value of epsilon
beta_c_0 = 2e-11;     % default value of chamber compressibility

% list of chamber aspect ratios. See benchmark case parameter sweep.
asp_list = [0.5; 1; 1.5; 2; 2.5; 3];

% list of chamber compressibilities for chambers height from 1 to 6 km. See
% benchmark case parameter sweep.
beta_c_list = [2.2526e-11; 2.5895e-11; 2.9388e-11;  3.3004e-11; 2e-11;  2.7243e-11];

% list of magma compressibilities. See benchmark case parameter sweep 
beta_m_list = [1e-10; 2.29264e-10; 3.65638e-10; 5.0949e-10; 5.77474e-10; 7.409320e-10];

beta_list = beta_m_list + beta_c_list;

% time
ts = linspace(0, 20, 1e3);

% epsilon, which can be determined from the regime diagram as a function of
% \omega R/c_s^r and \omega H/c_p^m
epsilon_list_1 = linspace(0, 2, 5);
epsilon_list_2 = [0, 1.5];

linetype_list = {'-', '--'}; % for plotting 

lambda_m = 1e10;      % magma Lame constant
lambda_r = 3e10;      % rock Lame constant
mu_m = 0;             % magma shear modulus
mu_r = 3e10;          % rock shear modulus
rho_m = 2.7e3;        % magma density
rho_r = 3e3;          % rock density
R_b = 1e3;            % caldera block radius
R_c = 1e3;            % cylindrical chamber radius
L_c = (2*R_c)*asp_0;  % cylindrical chamber height
L_b = 1e3;            % caldera block height
fs = 0.6;             % static friction coefficient
fd = 0.37;            % dynamic friction coefficient
sigma0 = 20e6;        % initial normal stress
tau0   = 12.01e6;     % initial shear stress
g = 9.8;              % gravitational constant
Dtau_str = -(fs-fd)*sigma0;    % fault strength drop

V_c_0 = pi*R_c^2*L_c;       % initial magma chamber volume
V_b = pi*R_b^2*L_b;         % caldera block volume
m_b = V_b * rho_r;          % caldera block mass
m_m = V_c_0 * rho_m;        % magma mass
m_prime = m_b + 1/3*m_m;    % effective inertial mass of combined caldera block and magma, assuming cylindrical chamber
K_m = lambda_m+ 2*mu_m/3;   % magma bulk modulus
K_r = lambda_r+ 2*mu_r/3;   % rock bulk modulus
beta_m = 1/K_m;             % magma compressibility
beta_0 = beta_m + beta_c_0; % default total compressibility (magma + chamber)

cp_m = sqrt((lambda_m + 2/3*mu_m)/rho_m); % P-wave velocity in magma

% the following is the analytical solution to Eqn. C5 (lumped model accounting for plane-wave radiation)
a_0 = m_prime;
b_0 = pi^2 * R_b^4/(beta_0*V_c_0);
d_0 = (2*pi*R_b*L_b)*Dtau_str;
c_0 = rho_m * cp_m * pi * R_c^2 * epsilon_0;

u = @(t, a, b, c, d) -d/(2*b*sqrt(-4*a*b + c^2)) * (-2*sqrt(-4*a*b+c^2) ...
    - c.*exp(0.5*(-c/a-sqrt(-4*a*b+c^2)/a).*t)...
    + sqrt(-4*a*b + c^2) .* exp(0.5*(-c/a-sqrt(-4*a*b+c^2)/a).*t)...
    + c.*exp(0.5*(-c/a+sqrt(-4*a*b+c^2)/a).*t)...
    + sqrt(-4*a*b + c^2) .* exp(0.5*(-c/a+sqrt(-4*a*b+c^2)/a).*t)...
    );

%% Varying degrees of radiation damping
cmap = autumn(length(asp_list));
figure;
for i = 1:length(epsilon_list_1)
    epsilon = epsilon_list_1(i);

    c = rho_m * cp_m * pi * R_c^2 * epsilon;
    s = real(u(ts, a_0, b_0, c, d_0)); % use real() to get spurious imaginary component

    % Numerically determine when slip ends
    v = gradient(s, ts);
    idc = intersect(find(abs(v) < 0.03), find(ts >0.5));

    if ~isempty(idc)
        idx = idc(1);
        tmax = ts(idx);
        s1 = s(ts <= tmax);
        s2 = s(ts > tmax);

        plot(ts(ts <= tmax), s1, 'Color', cmap(i, :), 'LineWidth', 2); hold on;
        plot(ts(ts > tmax), s2, 'Color', [cmap(i, :), 0.1], 'LineWidth', 2);
    else
        plot(ts, s, 'Color', cmap(i, :), 'LineWidth', 2); hold on;
    end
end
xlabel('Time (s)'); ylabel('Slip (m)')
xlim([0, 8]);

%% Varying chamber volume
cmap = copper(length(asp_list));

figure;
for i = 1:length(asp_list)
    asp = asp_list(i);

    L_c = (2*R_c)*asp;
    V_c = pi*R_c^2*L_c;
    m_m = V_c * rho_m; 
    m_prime = m_b + 1/3*m_m;

    for j = 1:length(epsilon_list_2)
        epsilon = epsilon_list_2(j);
  
        a = m_prime;
        b = pi^2 * R_b^4/(beta_0*V_c);
        c = rho_m * cp_m * pi * R_c^2 * epsilon;
        d = (2*pi*R_b*L_b)*Dtau_str;
        
        s = real(u(ts, a, b, c, d)); % get rid of small spurious imaginery components
        
        % Numerically determine when slip ends
        v = gradient(s, ts);
        idc = intersect(find(abs(v) < 0.03), find(ts >0.5));
        idx = idc(1);
        tmax = ts(idx);

        replace_idc = find(ts >= tmax);

        % replace the slip after tmax with flat line
        s(replace_idc) = repelem(s(replace_idc(1)), length(replace_idc));

        plot(ts, s, 'Color', cmap(i, :), 'LineWidth', 2, 'LineStyle', linetype_list{j}); hold on;
    end
end
xlabel('time (s)'); ylabel('Slip (m)')
ylim([-15, 0]);

%% Varying magma compressibility
cmap = copper(length(beta_m_list));

figure;
for i = 1:length(beta_m_list)
    beta_m = beta_m_list(i);
    beta = beta_m + beta_c_0;

    for j = 1:length(epsilon_list_2)
        epsilon = epsilon_list_2(j);
  
        b = pi^2 * R_b^4/(beta*V_c_0);
        c = rho_m * cp_m * pi * R_c^2 * epsilon;
        s = real(u(ts, a_0, b, c, d_0)); % get rid of small spurious imaginery components
        
        % determine when displacement end
        v = gradient(s, ts);
        idc = intersect(find(abs(v) < 0.03), find(ts >0.5));
        if ~isempty(idc)
            idx = idc(1);
            tmax = ts(idx);
            replace_idc = find(ts >= tmax);

            % replace the slip after tmax with flat line
            s(replace_idc) = repelem(s(replace_idc(1)), length(replace_idc));

        end

        plot(ts, s, 'Color', cmap(i, :), 'LineWidth', 2, 'LineStyle', linetype_list{j}); hold on;
    end
end
xlabel('time (s)'); ylabel('Slip (m)');
ylim([-15, 0]);

