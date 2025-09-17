% dispersion3D_FDFD.m
% This script computes the numerical dispersion characteristics of a 3D FDFD grid.
% It evaluates the relative error between the numerical wavenumber kn and analytical wavenumber k
% for a range of cells-per-wavelength (CPW) values, and plots error vs. incidence angle for a chosen CPW.

clear; clc; close all;

% Parameters
f = 1.5e9;
c = 3e8;
lambda = c/f;                % Wavelength (normalized)
k = 2*pi/lambda;             % Analytical wavenumber
CPW_list = [8, 10, 12, 16, 20, 24];  % Cells per wavelength to test
theta_range = linspace(0, pi/2, 181); % Incidence angles from 0 to 90 degrees
CPW_plot = 8;               % CPW for angle sweep plot

% Preallocate
error_CPW = zeros(size(CPW_list));

% Compute relative error for each CPW (isotropic propagation along x)
for idx = 1:length(CPW_list)
    CPW = CPW_list(idx);
    dx = lambda / CPW;
    % Numerical wavenumber along x-direction
    kn = (2/dx) * sin(k*dx/2);
    error_CPW(idx) = abs(kn - k) / k;
end

% Display results
fprintf('CPW    Relative error (kn vs k)\n');
for idx = 1:length(CPW_list)
    fprintf('%3d    %e\n', CPW_list(idx), error_CPW(idx));
end

% Plot error vs. CPW
figure;
semilogy(CPW_list, error_CPW, 'o-', 'LineWidth', 2);
xlabel('Cells per Wavelength (CPW)');
ylabel('Relative Error |kn - k| / k');
title('Dispersion Error vs. Grid Resolution');
grid on;

% Angle-dependent dispersion at chosen CPW_plot
dx_plot = lambda / CPW_plot;
error_theta = zeros(size(theta_range));
for t = 1:length(theta_range)
    theta = theta_range(t);
    % Analytical components
    kx = k * cos(theta);
    ky = k * sin(theta);
    kz = 0;
    % Numerical components using central differences
    knx = (2/dx_plot) * sin(kx*dx_plot/2);
    kny = (2/dx_plot) * sin(ky*dx_plot/2);
    knz = (2/dx_plot) * sin(kz*dx_plot/2);
    % Numerical wavenumber magnitude
    kn_mag = sqrt(knx^2 + kny^2 + knz^2);
    % Analytical magnitude
    k_mag = sqrt(kx^2 + ky^2 + kz^2);
    % Relative error
    error_theta(t) = (kn_mag ) / k_mag;
end

% Plot error vs. incidence angle
figure;
plot(theta_range*180/pi, error_theta, 'LineWidth', 2);
xlabel('Incidence Angle \theta (degrees)');
ylabel('Relative Error |kn - k| / k');
% title(sprintf('Dispersion Error vs. Angle (CPW = %d)', CPW_plot));
grid on;
