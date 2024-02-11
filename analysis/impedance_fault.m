%% Analyze impedance to 2D antiplane slip. 
% 
% The solution has applications to both ring fault slip and basaltic magma flow problem
% impedances are non-dimensionalized

%% Adjust font and line
set(0,'DefaultTextFontSize',16)
set(0,'DefaultAxesFontSize',16)

%% Define constants
N = 1e5;
gamma = double(eulergamma);
muN2muP = 1;
rhoN2rhoP = 1;
aR = linspace(1e-2, 1e2, N); % a = omega/c_n; use the dimensionless parameter on the negative side for plotting
%% Dimensionless impedance for mode k = 0
k = 0;

% compute impedance--------------------------------------------------------
z_n = @(aR, k) 1i.*(besselj(k-1,aR)./besselj(k,aR) - k./aR);
z_p = @(aR, k, muN2muP, rhoN2rhoP) -1i./muN2muP.*muN2muP.^(1/2).*(rhoN2rhoP).^(-1/2).*(besselh(k+1,aR.*muN2muP.^(1/2).*rhoN2rhoP.^(-1/2))./besselh(k,aR.*muN2muP.^(1/2).*(rhoN2rhoP).^(-1/2)) - k./aR.*muN2muP.^(-1/2).*(rhoN2rhoP).^(1/2));

Z0_n = z_n(aR, k);
Z0_p = z_p(aR, k, muN2muP, rhoN2rhoP);
Z0 = (Z0_p.*Z0_n)./(Z0_p+Z0_n);          % total impedance

% divide the impedance into real and imaginery parts--------------------
% Impedance inside of the ring fault
Z0_n_R1_idc = find(real(Z0_n) >= 0);
Z0_n_R2_idc = find(real(Z0_n) < 0);
aR_n_R1 = aR(Z0_n_R1_idc);
aR_n_R2 = aR(Z0_n_R2_idc);
Z0_n_R1 = real(Z0_n(Z0_n_R1_idc));% real positive
Z0_n_R2 = real(Z0_n(Z0_n_R2_idc));% real negative

Z0_n_I1_idc = find(imag(Z0_n) >= 0);
Z0_n_I2_idc = find(imag(Z0_n) < 0);
aR_n_I1 = aR(Z0_n_I1_idc);
aR_n_I2 = aR(Z0_n_I2_idc);
Z0_n_I1 = imag(Z0_n(Z0_n_I1_idc));% imaginery positive
Z0_n_I2 = imag(Z0_n(Z0_n_I2_idc));% imaginery negative

% Impedance outside of the ring fault
Z0_p_R1_idc = find(real(Z0_p) >= 0);
Z0_p_R2_idc = find(real(Z0_p) < 0);
aR_p_R1 = aR(Z0_p_R1_idc);
aR_p_R2 = aR(Z0_p_R2_idc);
Z0_p_R1 = real(Z0_p(Z0_p_R1_idc));% real positive
Z0_p_R2 = real(Z0_p(Z0_p_R2_idc));% real negative

Z0_p_I1_idc = find(imag(Z0_p) >= 0);
Z0_p_I2_idc = find(imag(Z0_p) < 0);
aR_p_I1 = aR(Z0_p_I1_idc);
aR_p_I2 = aR(Z0_p_I2_idc);
Z0_p_I1 = imag(Z0_p(Z0_p_I1_idc));% imaginery positive
Z0_p_I2 = imag(Z0_p(Z0_p_I2_idc));% imaginery negative

% total impedance
Z0_I1_idc = find(imag(Z0) >= 0);
Z0_I2_idc = find(imag(Z0) < 0);
Z0_R1_idc = find(real(Z0) >= 0);
Z0_R2_idc = find(real(Z0) < 0);
aR_I1 = aR(Z0_I1_idc);
aR_I2 = aR(Z0_I2_idc);
aR_R1 = aR(Z0_R1_idc);
aR_R2 = aR(Z0_R2_idc);
Z0_I1 = imag(Z0(Z0_I1_idc));% imaginery positive
Z0_I2 = imag(Z0(Z0_I2_idc));% imaginery negative
Z0_R1 = real(Z0(Z0_R1_idc));% real positive
Z0_R2 = real(Z0(Z0_R2_idc));% real negative

% absolute value of total impedance
Z0_abs = abs(Z0);
Z0_n_abs = abs(Z0_n);
Z0_p_abs = abs(Z0_p);

% compute limits-----------------------------------------------------------
% short wavelength limit
S0_n = repelem(nan, N);
S0_p = repelem(-muN2muP.^(-1).*muN2muP.^(1/2).*(rhoN2rhoP).^(-1/2), N);

% long wavelength limit
L0_n = -1i.*aR/2;
L0_p = 2i./muN2muP./aR./(2*gamma-1i.*pi+2.*log(aR./2.*muN2muP.^(1/2).*rhoN2rhoP.^(-1/2)));

%% Impedance amplitude for k = 0 mode
figure;
loglog(aR, Z0_n_abs, 'k--', 'LineWidth', 2); hold on;
loglog(aR, Z0_p_abs, 'k:', 'LineWidth', 2);
loglog(aR, Z0_abs, 'k-', 'LineWidth', 2);
loglog(aR(aR >= 1), abs(S0_n(aR >= 1)), 'b--', 'LineWidth',2);
loglog(aR(aR >= 1), abs(S0_p(aR >= 1)), 'b:',  'LineWidth',2);
loglog(aR(aR <= 1), abs(L0_n(aR <= 1)), 'r--', 'LineWidth', 2);
loglog(aR(aR <= 1), abs(L0_p(aR <= 1)), 'r:', 'LineWidth', 2);
text1 = text(0.15, 10^(-1), '$|\frac{-i \omega \rho R}{2}|$', 'Interpreter','latex', 'FontSize', 16);
text2 = text(0.06, 10^(1), '$|\frac{i \mu_{+}}{c_s^{+} (\omega R/c_s^{+}) \log(\omega R/(2 c_s^{+})) + \omega R (\gamma - i\pi/2)}|$', 'Interpreter','latex', 'FontSize', 16);
text3 = text(12, 10^(1), '$|\frac{\mu_{+}}{c_s^{+}}|$', 'Interpreter','latex', 'FontSize', 16);
set(text1, 'Color', 'r');
set(text2, 'Color', 'r');
set(text3, 'Color', 'b');
legend('$|\tilde{Z}^{-}|$', '$|\tilde{Z}^{+}|$', '$|\tilde{Z}|$','Interpreter','latex');
xlabel('\omega R/c_s^{-}'); 
grid on;
title(strcat('$k = 0, \mu^{-}/\mu^{+} =',num2str(muN2muP), ', \rho^{-}/\rho^{+} = 1$'), 'Interpreter', 'latex')

%% Real and imaginery components of the impedance for k = 0 mode
% Notes:
% Z0 does not have a real positive component

figure;
subplot(2, 1, 1)
%h1 = loglog(aR_R1, Z0_R1 , 'k-', 'LineWidth', 2); hold on;                  % real positive (does not exist)
h2 = loglog(aR_I1, Z0_I1 , '-',  'LineWidth', 2, 'Color', [0.6 0.6 0.6]);    % imaginery positive
xlim([10^(-2), 10^2]); ylim([10^-3, 10^1]);
grid on;
subplot(2, 1, 2)
h3 = loglog(aR_R2, Z0_R2 , 'k-', 'LineWidth', 2); hold on;                   % real negative
h4 = loglog(aR_I2, Z0_I2 , '-', 'LineWidth', 2, 'Color', [0.6 0.6 0.6]);     % imaginery negative
loglog(aR(aR >= 1), S0_p(aR >= 1), 'b-', 'LineWidth', 2);
loglog(aR(aR < 1),  imag(L0_n(aR < 1)), 'r-', 'LineWidth', 2);
ylim([-10^6, -10^(-3)]); 
legend([h3, h4], 'Re[$\tilde{Z}_{0}$]', 'Im[$\tilde{Z}_{0}$]','Interpreter', 'latex')
text1 = text(10^(1), -10^2, '$-\frac{\mu_{+}}{c_s^{+}}$', 'Interpreter','latex', 'FontSize', 16);
text2 = text(10^(-1), -10^2, 'Im[$(\frac{-i \omega \rho R}{2})]$', 'Interpreter','latex', 'FontSize', 16);
set(text1, 'Color', 'b');
set(text2, 'Color', 'r');
xlabel('log_{10} \omega R/c_s^{-}'); 
sgtitle('Axisymmetric mode')
grid on;

%{
% Real component of k = 0 mode

% Notes:
% Z0_n does not have (negative or positive) real components (therefore purely energy storage)
% Z0_p does not have positive real components 
% Z0 does not have positive real components

% Black lines for positive component
% Gray lines for negative component

figure;
h1 = loglog(aR_n_R1, Z0_n_R1 , 'k--', 'LineWidth', 2); hold on;
h2 = loglog(aR_n_R2, -Z0_n_R2 , '--', 'LineWidth', 2, 'Color', [0.6 0.6 0.6]);
h3 = loglog(aR_p_R1, Z0_p_R1, 'k:', 'LineWidth', 2);
h4 = loglog(aR_p_R2, -Z0_p_R2 , ':', 'LineWidth', 2, 'Color', [0.6 0.6 0.6]);
h5 = loglog(aR_R1, Z0_R1, 'k-', 'LineWidth', 2);
h6 = loglog(aR_R2, -Z0_R2 , '-', 'LineWidth', 2, 'Color', [0.6 0.6 0.6]);
loglog(aR(aR >= 1), -real(S0_p(aR >= 1)), 'b:',  'LineWidth',2);
%legend([h1, h2, h3, h4, h5, h6], '$Re(\tilde{Z}^{-}) >= 0$', '$Re(\tilde{Z}^{-}) < 0$', '$Re(\tilde{Z}^{+}) >= 0$', '$Re(\tilde{Z}^{+}) < 0$', '$Re(\tilde{Z}) >= 0$', '$Re(\tilde{Z}) < 0$', 'Interpreter','latex');
legend([h4, h6], 'Re$[\tilde{Z}^{+}$] $< 0$', 'Re[$\tilde{Z}^{-}$] $ < 0$', 'Interpreter','latex');
xlabel('\omega R/c_s^{-}'); ylabel('log_{10} Pa s m^{-1}')
grid on;
title('k = 0, real component')

% Imaginary component of k = 0 mode

% Notes:
% Z0_p does not have positive imaginery components (therefore does not have energy storage)

% Black lines for positive component
% Gray lines for negative component

figure;
h1 = loglog(aR_n_I1, Z0_n_I1 , 'k--', 'LineWidth', 2); hold on;
h2 = loglog(aR_n_I2, -Z0_n_I2 , '--', 'LineWidth', 2, 'Color', [0.6 0.6 0.6]);
h3 = loglog(aR_p_I1, Z0_p_I1, 'k:', 'LineWidth', 2);
h4 = loglog(aR_p_I2, -Z0_p_I2 , ':', 'LineWidth', 2, 'Color', [0.6 0.6 0.6]);
h5 = loglog(aR_I1, Z0_I1, 'k-', 'LineWidth', 2);
h6 = loglog(aR_I2, -Z0_I2 , '-', 'LineWidth', 2, 'Color', [0.6 0.6 0.6]);
loglog(aR(aR <= 1), -imag(L0_n(aR <= 1)), 'r--', 'LineWidth', 2);
loglog(aR(aR <= 1), -imag(L0_p(aR <= 1)), 'r:', 'LineWidth', 2);
%legend([h1, h2, h3, h4, h5, h6], '$Im(\tilde{Z}^{-}) >= 0$', '$Im(\tilde{Z}^{-}) < 0$', '$Im(\tilde{Z}^{+}) >= 0$', '$Im(\tilde{Z}^{+}) < 0$', '$Im(\tilde{Z}) >= 0$', '$Im(\tilde{Z}) < 0$', 'Interpreter','latex');
legend([h1, h2, h4, h5, h6], 'Im[$\tilde{Z}^{-}]$ $>= 0$', 'Im[$\tilde{Z}^{-}$] $< 0$', 'Im[$\tilde{Z}^{+}$] $< 0$', 'Im$[\tilde{Z}]$ $>= 0$', 'Im$[\tilde{Z}]$ $< 0$', 'Interpreter','latex');
xlabel('\omega R/c_s^{-}'); ylabel('log_{10} Pa s m^{-1}')
grid on;
title('k = 0, imaginary component')
%}


%% Plot the slip profile for k = 0, k = 1, and the asymmetric mode
theta = linspace(0, 2*pi, N);

figure;
plot(theta./pi, -repelem(1, N), 'k-', 'LineWidth', 1); hold on;
plot(theta./pi, -cos(theta), '-', 'LineWidth', 1, 'Color', [0.6 0.6 0.6]);
plot(theta./pi, 0.5*(-cos(theta)-1), '-', 'LineWidth', 1);
text(0.7, 1.1, 'Re [$\hat{\delta}_{1} (\theta, \omega)/D(\omega)$]','Interpreter','latex')
text(0.4, 0.2, 'Re [$\hat{\delta}_{asymmetric}(\theta, \omega)/D(\omega) = \frac{1}{2} (\hat{\delta}_0 + \hat{\delta}_1)/D(\omega)]$','Interpreter','latex')
text(0.7, -0.8, 'Re [$\hat{\delta}_{0}(\theta, \omega)/D(\omega))$]','Interpreter','latex')
ylim([-1.25, 1.25]);xlim([0, 2]);

%% Dimensionless impedance for k = 1
% Notes:
% 1. The amplitude of total impedance at the low frequency limit is higher
% than the harmonic sum of the asymptotes inside and outside of the ring
% fault, because the asymptotes are PURELY imaginery, but the actual
% impedances are complex. 

k = 1;
muN2muP_k1 = 1;

Z1_n = z_n(aR, k);
Z1_p = z_p(aR, k, muN2muP_k1, rhoN2rhoP);
Z1 = (Z1_p.*Z1_n)./(Z1_p+Z1_n);  % total impedance

Z1_n_abs = abs(Z1_n);
Z1_p_abs = abs(Z1_p);
Z1_abs = abs(Z1);

L1_n = 1i./aR;
L1_p = 1./aR./muN2muP_k1./1i;
S1_p = repelem(-muN2muP_k1.^(-1).*muN2muP_k1.^(1/2).*(rhoN2rhoP).^(-1/2), N);

figure;
loglog(aR, Z1_n_abs, 'k--', 'LineWidth', 2); hold on;
loglog(aR, Z1_p_abs, 'k:', 'LineWidth', 2);
loglog(aR, Z1_abs, 'k-', 'LineWidth', 2);
loglog(aR(aR <= 1), abs(L1_n(aR <= 1)), 'r--', 'LineWidth', 2);
loglog(aR(aR <= 1), abs(L1_p(aR <= 1)), 'r:', 'LineWidth', 2);
loglog(aR(aR >= 1), abs(S1_p(aR >= 1)), 'b:', 'LineWidth', 2);
text1 = text(10^(-1), 10^0, '$|-\frac{\mu_{-}}{i \omega R}|$', 'Interpreter','latex', 'FontSize', 16);
text2 = text(10^(-1), 10^2, '$|\frac{\mu_{+}}{i \omega R}|$', 'Interpreter','latex', 'FontSize', 16);
text3 = text(10^(1), 10^1, '$|-\frac{\mu_{+}}{c_s^{+}}|$', 'Interpreter','latex', 'FontSize', 16);
set(text1, 'Color', 'r'); set(text2, 'Color', 'r'); set(text3, 'Color', 'b');
title(strcat('$k = 1, \mu^{-}/\mu^{+} =',num2str(muN2muP_k1), ', \rho^{-}/\rho^{+} = 1$'), 'Interpreter', 'latex')
xlabel('\omega R/c_s^{-}');
grid on;

%% Asymmetric mode
% Notes:
% 1. The imaginery component of total impedance dominates over the real component in the low frequency limit.
% 2. The real component of total impedance dominates over the imaginery component in the low frequency limit.


% composite dimensionless impedance for asymmetric collapse
Zc = (Z0+Z1)./2; 

% divide the impedance into real and imaginery parts--------------------

% total impedance
Zc_I1_idc = find(imag(Zc) >= 0);
Zc_I2_idc = find(imag(Zc) < 0);
Zc_R1_idc = find(real(Zc) >= 0);
Zc_R2_idc = find(real(Zc) < 0);
aR_I1 = aR(Zc_I1_idc);
aR_I2 = aR(Zc_I2_idc);
aR_R1 = aR(Zc_R1_idc);
aR_R2 = aR(Zc_R2_idc);
Zc_I1 = imag(Zc(Zc_I1_idc));% imaginery positive
Zc_I2 = imag(Zc(Zc_I2_idc));% imaginery negative
Zc_R1 = real(Zc(Zc_R1_idc));% real positive
Zc_R2 = real(Zc(Zc_R2_idc));% real negative

% 1. Impedance amplitude for asymmetric mode, compared to those of axisymmetric and tilting mode
figure;
loglog(aR, abs(Z0), 'k-', 'LineWidth', 2); hold on;
loglog(aR, abs(Z1), '-', 'LineWidth', 2,  'Color', [0.6 0.6 0.6]);
loglog(aR, abs(Zc), 'LineWidth', 2);
ylim([10^-6, 10^6])
xlabel('log_{10} \omega R/c_s^{-}'); 
grid on;
legend('$|\tilde{Z}_0|$', '$|\tilde{Z}_1|$', '$|\tilde{Z}_{asymmetric}|$','Interpreter','latex');
%}

% Asymptotic limits for total impedance 
Sa = S1_p;
La = 2.*1i./(1+4*gamma-2.*1i*pi - 4.*log(2) + 4*log(aR))./aR.^3;
La_real = real(La);
La_imag = imag(La);

% 2. Impedance in terms of real and imaginery components
figure;
subplot(2, 1, 1)
%h1 = loglog(aR_R1, Zc_R1 , 'k-', 'LineWidth', 2); hold on;    % real positive (does not exist)
h2 = loglog(aR_I1, Zc_I1 , '-',  'LineWidth', 2, 'Color', [0.6 0.6 0.6]);                                    % imaginery positive
xlim([10^(-2), 10^2]); ylim([10^-3, 10^1]);
grid on;
subplot(2, 1, 2)
h3 = loglog(aR_R2, Zc_R2 , 'k-', 'LineWidth', 2); hold on;      % real negative
h4 = loglog(aR_I2, Zc_I2 , '-', 'LineWidth', 2, 'Color', [0.6 0.6 0.6]);                                     % imaginery negative
loglog(aR(aR >= 1), Sa(aR >= 1)./2, 'b-', 'LineWidth', 2);
loglog(aR(aR < 1),  La_real(aR < 1), 'r--', 'LineWidth', 2);
loglog(aR(aR < 1),  La_imag(aR < 1), 'r--', 'LineWidth', 2);
ylim([-10^6, -10^(-3)]); 
legend([h3, h4], 'Re$[\tilde{Z}_{asymmetric}]$', 'Im$[\tilde{Z}_{asymmetric}]$','Interpreter', 'latex')
text1 = text(10^(1), -10^2, '$-\frac{\mu_{+}}{2 c_s^{+}}$', 'Interpreter','latex', 'FontSize', 16);
text2 = text(10^(-1), -10^2, '$\frac{2 i \mu_{-}}{c_s^{-}(1 + 4 \gamma - 2 \pi i + 4 \log \frac{\omega R}{2 c_s^{-}})(\frac{\omega R}{c_s^{-}})^{3})}$', 'Interpreter','latex', 'FontSize', 16);
set(text1, 'Color', 'b');
set(text2, 'Color', 'r');
xlabel('log_{10} \omega R/c_s^{-}'); 
sgtitle('Asymmetric mode')
grid on;



%% Impedance for various shear moduli contrasts (k = 0 mode)
muN2muP_list = [1e-2, 1e-1, 1, 1e1, 1e2];
rhoN2rhoP_list = [1e-2, 1e-1, 1, 1e1, 1e2];

figure;
for i = 1:length(muN2muP_list)
    
    Z0_n = z_n(aR, 0);
    Z0_p = z_p(aR, 0, muN2muP_list(i), rhoN2rhoP);
    Z0 = (Z0_p.*Z0_n)./(Z0_p+Z0_n);

    % absolute value of impedance
    Z0_abs = abs(Z0);
    Z0_n_abs = abs(Z0_n);
    Z0_p_abs = abs(Z0_p);

    loglog(aR, Z0_abs, 'k-', 'LineWidth', 2); hold on;
end
xlim([2, 6]); ylim([10^(-5), 10^2]);
xlabel('\omega R/c_s^{-}'); 
title('Effect of shear modulus contrasts on k = 0 mode')

figure;
for i = 1:length(rhoN2rhoP_list)

    Z0_n = z_n(aR, 0);
    Z0_p = z_p(aR, 0, muN2muP, rhoN2rhoP_list(i));
    Z0 = (Z0_p.*Z0_n)./(Z0_p+Z0_n);

    % absolute value of impedance
    Z0_abs = abs(Z0);
    Z0_n_abs = abs(Z0_n);
    Z0_p_abs = abs(Z0_p);

    loglog(aR, Z0_abs, 'k-', 'LineWidth', 2); hold on;
end
xlim([2, 6]); ylim([10^(-5), 10^2]);
xlabel('\omega R/c_s^{-}'); 
title('Effect of density contrasts on k = 0 mode')






