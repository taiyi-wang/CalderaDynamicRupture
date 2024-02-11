%% Inversion procedure to fit the desired relaxation functions for Maxwell rheology

current_dir  = pwd;
idcs   = strfind(current_dir ,'/');
above_dir = current_dir(1:idcs(end)-1); % one level above current directory

addpath(fullfile(above_dir, 'src_code/.'))

%% Set up target relaxation functions
Nsteps = 3e3;    % number of iterations for the inversion

% choose fixed parameters--------------------------------
n = 9;           % number of mechanisms
vp_0 = 1e3;      % p wave velocity at infinite frequency (place holder values)
vs_0 = 1e3;      % s wave velocity at infinite frequency (place holder values)

% desired magma material properties-----------------------------
% best fit so far (visually)
% 1. mu = 1e7, eta = 1e8, rho = 2700, K = 1e10 (silicic magma; default values)
% 2. mu = 1e5, eta = 1e5; rho = 2200, K = 1e9; (intermediate vesicular magma; default values)

mu = 1e7;         % shear modulus (See Okmura et al., 2010; for silicate melt without bubbles)
eta = 1e8;        % viscosity (Okmura et al., 2010: 0.4-0.5 wt% of water at 830 degree C)
rho = 2700;       % density
K = 1e10;         % bulk modulus (Huppert and Woods 2002, as cited in Gottsmann, 2020)

% compute desired relaxation functions
G1_M = @(t) 2*mu*exp(-t/(eta/mu));
G2_M = 3*K;

N = 1e4;
tau_M = eta/mu;
t = linspace(0, 20, N)'.*tau_M;

figure;
yyaxis left
h1 = plot(t./tau_M , G1_M(t)./1e9, 'r--', 'LineWidth', 2);
ylabel('GPa')
yyaxis right
h2 = yline(G2_M./1e9, 'b--', 'LineWidth', 2);
grid on;
legend([h1, h2], 'G_1^{M}', 'G_2^{M}')
xlabel('t/\tau^{M}'); ylabel('GPa')
ax = gca;
ax.YAxis(1).Color = 'r';
ax.YAxis(2).Color = 'b';

%% Inversion
% setup range of inversion parameters----------------------------------------

% parameters to vary
% in the order of:
% log10(f_c), log10(f_ratio), log10(Qp), log10(Qs), log10(vp_t), log10(vs_t)
var_bnds = [-1, 1; 2, 4; 4, 6; 0, 2; 2, 3.5; 0, 2.5];
params = [n, rho, eta, mu, K, vp_0, vs_0];

rng default
f1 = @(md) MaxwellObjective(md, params, t); 
lb = var_bnds(:,1);ub = var_bnds(:,2);
options = optimoptions('surrogateopt','MaxFunctionEvaluations',Nsteps, 'PlotFcn','surrogateoptplot');
[md,fval] = surrogateopt(f1,lb,ub, options);


%% Plot best fit model
[G1, G2, params] = cmp_relaxationFunc(10.^md(3), 10.^md(4), n, 10.^md(1), 10.^md(2), vp_0, vs_0, 10.^md(5), 10.^md(6), rho);

figure;
yyaxis left
h1 = loglog(t./tau_M , G1(t')', 'r-', 'LineWidth', 2); hold on;
h2 = loglog(t./tau_M , G1_M(t')', 'r--', 'LineWidth', 2);
ylim([10^(-2), 10^(10)]);
ylabel('log_{10} [Pa]')
yyaxis right
h3 = loglog(t./tau_M , G2(t')', 'b-','LineWidth', 2); 
h4 = yline(G2_M, 'b--', 'LineWidth', 2);
%ylim([2.5e10, 3.5e10])
ylim([2.5e9, 3.5e9])
legend([h1, h2, h3, h4], 'G_1', 'G_1^{M}', 'G_2', 'G_2^{M}')
xlabel('log_{10} t/\tau^{M}'); ylabel('log_{10} [Pa]')
grid on;
ax = gca;
ax.YAxis(1).Color = 'r';
ax.YAxis(2).Color = 'b';
ax.GridColor = [1, 1, 1];
