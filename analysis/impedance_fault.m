%% Analyze impedance to antiplane slip in Fourier domain. 

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
Z0_R1 = real(Z0(Z0_R1_idc));% imaginery positive
Z0_R2 = real(Z0(Z0_R2_idc));% imaginery negative

% absolute value of impedance
Z0_abs = abs(Z0);
Z0_n_abs = abs(Z0_n);
Z0_p_abs = abs(Z0_p);

% compute limits-----------------------------------------------------------
% short wavelength limit
S0_n = repelem(nan, N);
S0_p = repelem(-muN2muP.^(-1).*muN2muP.^(1/2).*(rhoN2rhoP).^(-1/2), N);

% long wavelength limit
L0_n = -1i.*aR/2;
L0_p = 2i./muN2muP.*muN2muP.^(1/2).*rhoN2rhoP.^(-1/2)./(2*gamma-1i.*pi+2.*log(aR./2.*muN2muP.^(1/2).*rhoN2rhoP.^(-1/2)))./(muN2muP.^(1/2).*rhoN2rhoP.^(-1/2))./aR;


%% Absolute Impedance for K = 0 mode
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
legend('$|\hat{Z}^{-}|$', '$|\hat{Z}^{+}|$', '$|\hat{Z}|$','Interpreter','latex');
xlabel('\omega R/c_s'); 
grid on;
title('k = 0')

%% Real component of k = 0 mode
figure;
h1 = loglog(aR_n_R1, Z0_n_R1 , 'k--', 'LineWidth', 2); hold on;
loglog(aR_n_R2, -Z0_n_R2 , '--', 'LineWidth', 2, 'Color', [0.6 0.6 0.6]);
h2 = loglog(aR_p_R1, Z0_p_R1, 'k:', 'LineWidth', 2);
loglog(aR_p_R2, -Z0_p_R2 , ':', 'LineWidth', 2, 'Color', [0.6 0.6 0.6]);
h3 = loglog(aR_R1, Z0_R1, 'k-', 'LineWidth', 2);
loglog(aR_R2, -Z0_R2 , '-', 'LineWidth', 2, 'Color', [0.6 0.6 0.6]);
loglog(aR(aR >= 1), -real(S0_p(aR >= 1)), 'b:',  'LineWidth',2);
text1 = text(12, 10^7, '$-Re(\frac{\mu_{+}}{c_s^{+}})$', 'Interpreter','latex', 'FontSize', 16);
set(text1, 'Color', 'b');
legend([h1, h2, h3], '$Re(\hat{Z}^{-})$', '$-Re(\hat{Z}^{+})$', '$Re(\hat{Z})$','Interpreter','latex');
xlabel('\omega R/c_s'); ylabel('log_{10} Pa s m^{-1}')
grid on;
title('k = 0, real component')

%% Imaginary component of k = 0 mode
figure;
h1 = loglog(aR_n_I1, Z0_n_I1 , 'k--', 'LineWidth', 2); hold on;
loglog(aR_n_I2, -Z0_n_I2 , '--', 'LineWidth', 2, 'Color', [0.6 0.6 0.6]);
h2 = loglog(aR_p_I1, Z0_p_I1, 'k:', 'LineWidth', 2);
loglog(aR_p_I2, -Z0_p_I2 , ':', 'LineWidth', 2, 'Color', [0.6 0.6 0.6]);
h3 = loglog(aR_I1, Z0_I1, 'k-', 'LineWidth', 2);
loglog(aR_I2, -Z0_I2 , '-', 'LineWidth', 2, 'Color', [0.6 0.6 0.6]);
loglog(aR(aR <= 1), -imag(L0_n(aR <= 1)), 'r--', 'LineWidth', 2);
loglog(aR(aR <= 1), -imag(L0_p(aR <= 1)), 'r:', 'LineWidth', 2);
text1 = text(0.15, 10^5, '$-Im(\frac{-i \omega \rho R}{2})$', 'Interpreter','latex', 'FontSize', 16);
text2 = text(0.06, 10^8, '$-Im(\frac{i \mu_{+}}{c_s^{+} (\omega R/c_s^{+}) \log(\omega R/(2 c_s^{+})) + \omega R (\gamma - i\pi/2)})$', 'Interpreter','latex', 'FontSize', 16);
set(text1, 'Color', 'r'); set(text2, 'Color', 'r');
legend([h1, h2, h3], '$Im(\hat{Z}^{-})$', '$-Im(\hat{Z}^{+})$', '$Im(\hat{Z})$','Interpreter','latex');
xlabel('\omega R/c_s'); ylabel('log_{10} Pa s m^{-1}')
grid on;
title('k = 0, imaginary component')


%% Plot the modes
theta = linspace(0, 2*pi, N);

figure;
plot(theta./pi, -repelem(1, N), 'k-', 'LineWidth', 1); hold on;
plot(theta./pi, -cos(theta), '-', 'LineWidth', 1, 'Color', [0.6 0.6 0.6]);
plot(theta./pi, 0.5*(-cos(theta)-1), '-', 'LineWidth', 1);
text(1, 1.1, '$|i \omega \hat{\delta}_{1}|$','Interpreter','latex')
text(1, 0.5, '$|i \omega \hat{\delta}_{asymmetric}|$','Interpreter','latex')
text(1, -0.8, '$|i \omega \hat{\delta}_{0}|$','Interpreter','latex')
ylim([-1.25, 1.25]);xlim([0, 2]);


%% Impedance for mixed, asymmetric mode
k = 1;

Z1_n = z_n(aR, k);
Z1_p = z_p(aR, k, muN2muP, rhoN2rhoP);
Z1 = (Z1_p.*Z1_n)./(Z1_p+Z1_n);  % total impedance

Zc = 2.*Z0.*Z1./(Z0+Z1); % composite impedance for asymmetric collapse

figure;
loglog(aR, abs(Z0), 'k-', 'LineWidth', 2); hold on;
loglog(aR, abs(Z1), '-', 'LineWidth', 2,  'Color', [0.6 0.6 0.6]);
loglog(aR, abs(Zc), 'LineWidth', 2);
loglog(aR, abs(L0_n), 'r--', 'LineWidth', 2);
loglog(aR, abs(S0_p), 'b--', 'LineWidth', 2);
ylim([10^-6, 10^6])
xlabel('log_{10} \omega R/c_s'); ylabel('log_{10} Pa s m^{-1}');
grid on;
legend('$|\hat{Z}_0|$', '$|\hat{Z}_1|$', '$|\hat{Z}_{asymmetric}|$','Interpreter','latex');


%% Impedance for various shear moduli contrasts
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
xlabel('\omega R/c_s'); 
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
xlabel('\omega R/c_s'); 
title('Effect of density contrasts on k = 0 mode')






