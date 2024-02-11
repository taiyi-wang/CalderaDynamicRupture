%% Script to demonstrate the effects of neglecting radiation damping in rate-and-state based lumped model
% The model and code is adopted from Segall & Anderson, 2021 (PNAS)

current_dir  = pwd;
idcs   = strfind(current_dir ,'/');
above_dir = current_dir(1:idcs(end)-1); % one level above current directory

addpath(fullfile(above_dir, 'src_code/.'))
addpath(fullfile(above_dir, 'data/.'))

% Load sample model
load SampleMod.mat

%% Run
% 1. without radiation damping-------------------------------------------
[dhat, cfail, t, m, const] = predicted_data(mlm, const, 0);

% Dimensionalize results
% define constants
d_c = const.d_c; A = const.A; B = const.B; sigma = const.sigma; f0 = const.f0;
kappa = const.kappa; alpha = const.alpha; eta = const.eta; pout = const.pout;  
tstar = const.tstar; L = const.L; pstar = const.pstar; 

tau = 1e-6*A*sigma*asinh(0.5*m(:,3).*exp(f0/A) ...
       .*(m(:,2)/d_c).^(B/A));

t1 = t.*tstar.*86400;              % time (s)
p1 = m(:,1)*pstar;                 % chamber-averaged  pressure (Pa)
v1= m(:,3)*L/(tstar*86400);        % slip rate (m/s)
d1 = -m(:,4)*L;                    % slip (m)
tau1 = real(tau)*pstar.*1e6;       % ring fault-averaged shear stress (Pa)

% 2. with radiation damping----------------------------------------------
[dhat, cfail, t, m, const] = predicted_data(mlm, const);

% Dimensionalize results
% define constants
d_c = const.d_c; A = const.A; B = const.B; sigma = const.sigma; f0 = const.f0;
kappa = const.kappa; alpha = const.alpha; eta = const.eta; pout = const.pout;  
tstar = const.tstar; L = const.L; pstar = const.pstar; 

tau = 1e-6*A*sigma*asinh(0.5*m(:,3).*exp(f0/A) ...
       .*(m(:,2)/d_c).^(B/A));

t2 = t.*tstar.*86400;              % time (s)
p2 = m(:,1)*pstar;                 % chamber-averaged  pressure (Pa)
v2= m(:,3)*L/(tstar*86400);        % slip rate (m/s)
d2 = -m(:,4)*L;                    % slip (m)
tau2 = real(tau)*pstar.*1e6;       % ring fault-averaged shear stress (Pa)


figure;
subplot(2, 2, 1)
plot(t1./86400, p1./1e6, 'r-', 'LineWidth', 1); hold on;
plot(t2./86400, p2./1e6, 'k-', 'LineWidth', 1);
subplot(2, 2, 2)
semilogy(t1./86400, v1, 'r-', 'LineWidth', 1); hold on;
semilogy(t2./86400, v2, 'k-', 'LineWidth', 1);
subplot(2, 2, 3)
plot(t1./86400, tau1./1e6, 'r-', 'LineWidth', 1);  hold on;
plot(t2./86400, tau2./1e6, 'k-', 'LineWidth', 1);
subplot(2, 2, 4)
plot(t1./86400, d1, 'r-', 'LineWidth', 1);  hold on;
plot(t2./86400, d2, 'k-', 'LineWidth', 1);

%% phase diagram
v0 = L/(tstar*86400);
v_list = 10.^linspace(-6, 0, 100);
f_ss = @(v) f0 + (A-B) .* log(v./v0);
signom = pstar*sigma;

figure;
semilogx(v1, tau1./signom , 'r-', 'LineWidth', 1); hold on;
semilogx(v2, tau2./signom , 'k-', 'LineWidth', 1);
semilogx(v1(1), tau1(1)./signom , 'ro', 'MarkerSize', 10);
semilogx(v1(1), tau2(1)./signom , 'ko', 'MarkerSize', 10);
semilogx(v_list, f_ss(v_list) , '-', 'Color', [0, 0, 0, 0.5], 'LineWidth', 1);
grid on;
xlabel('log_{10} v'); ylabel('\tau/\sigma'); 




