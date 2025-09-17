% plot_dispersion.m
clear; clc;

%% === user parameters ===
f = 1.5e9;
c = 3e8;
lambda = c/f;                % Wavelength (normalized)
k = 2*pi/lambda;             % Analytical wavenumber
CPW_list = [5, 8, 10, 12, 14,16,18, 20, 24,28, 30,35, 40,45, 50];
M_list = [1, 2, 3, 4];
CPW_plot_theta = 5;
%% error vs. incidence angle
% theta = [1:0.5:90]*pi/180;
% for t = 1: length(theta)
%     k_n_1 = dispersion(lambda, M_list(1), CPW_list(4), theta(t), k);
%     k_n_2 = dispersion(lambda, M_list(2), CPW_plot_theta, theta(t), k);
%     error_1(t) = abs(k_n_1-k)/k;
%     error_2(t) = abs(k_n_2-k)/k;
% 
% end
% 
% 
% % Plot error vs. incidence angle
% figure;
% plot(theta, error_1, 'LineWidth', 2);hold on;
% plot(theta, error_2, 'LineWidth', 2);hold on;
% xlabel('Incidence Angle \theta',FontSize=12);
% ylabel('Relative Error |kn-k|/k',FontSize=12);
% legend("M=1(FDFD),CPW= 15","M=2,CPW=5",Location="north",Fontsize=12);
% title(sprintf('Dispersion Error vs. Angle'),FontSize=12);
% grid on;

%% error vs CPW
for idx = 1:length(CPW_list)
    k_n_1 = dispersion(lambda, M_list(1), CPW_list(idx), 45*pi/180, k);
    k_n_2 = dispersion(lambda, M_list(2), CPW_list(idx), 45*pi/180, k);
    k_n_3 = dispersion(lambda, M_list(3), CPW_list(idx), 45*pi/180, k);
    k_n_4 = dispersion(lambda, M_list(4), CPW_list(idx), 45*pi/180, k);

    error_1(idx) = abs(k_n_1-k)/k;
    error_2(idx) = abs(k_n_2-k)/k;
    error_3(idx) = abs(k_n_3-k)/k;
    error_4(idx) = abs(k_n_4-k)/k;
end

% Plot error vs CPW
figure;
semilogy(CPW_list, error_1, 'o-', 'LineWidth', 2);hold on;
semilogy(CPW_list, error_2, 'o-', 'LineWidth', 2);hold on;
semilogy(CPW_list, error_3, 'o-', 'LineWidth', 2);hold on;
semilogy(CPW_list, error_4, 'o-', 'LineWidth', 2);hold on;
xlabel('Cells per Wavelength (CPW)',FontSize=12);
ylabel('Relative Error |kn - k| / k',FontSize=12);
legend("M=1(FDFD)","M=2","M=3","M=4","Location","best",Fontsize=12);
title('Dispersion Error vs. Grid Resolution',FontSize=12);
grid on;

%% === numerical dispersion of DSCFD, 2-D case ====
function k_n = dispersion(lambda ,M, CPW, theta, k)
% === dsc parameters ===
if M == 1
    coeffs = [1];
    coor = [0.5];
elseif M == 2
    coeffs = [9/8, -1/24];
    coor = [0.5, 1.5];
elseif M == 3
    coeffs = [75/64, -25/384, 3/640];
    coor = [0.5, 1.5, 2.5];
elseif M == 4
    coeffs = [1225/1024, -245/3072, 49/5120, -5/7168]; % x-4 to x+4
    coor = [0.5, 1.5, 2.5, 3.5]; 
end

delta_x = lambda/CPW;
delta_y = lambda/CPW; 
k_x = k*cos(theta); k_y = k*sin(theta);
knx = (2/delta_x)*sum(coeffs.*sin(k_x*coor*delta_x));
kny = (2/delta_y)*sum(coeffs.*sin(k_y*coor*delta_y));
k_n = (knx^2+kny^2)^0.5;
end