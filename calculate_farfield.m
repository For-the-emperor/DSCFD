%CALCULATE_FARFIELD 用于计算远场RCS
% 定义传输边界

li = nc_farbuffer+1;
lj = nc_farbuffer+1;
lk = nc_farbuffer+1;
ui = Nx-nc_farbuffer+1;
uj = Ny-nc_farbuffer+1;
uk = Nz-nc_farbuffer+1;
% 计算外推边界的电流J和磁流M
calculate_J_and_M;
%% 计算和显示远场
disp("postprocessing and displaying simulation result");
calculate_radiated_power;
j = sqrt(-1);
number_of_angles = 180;

% parameters used by polar plotting functions
step_size = 10;       % increment between the rings in the polar grid
Nrings = 4;           % number of rings in the polar grid
line_style1 = 'b-';   % line style for theta component
line_style2 = 'r--';  % line style for phi component
scale_type = 'dB';    % linear or dB

plot_type = 'RCS';
%% xy plane
% % ===============================================
% farfield_theta = zeros(number_of_angles, 1);
% farfield_phi   = zeros(number_of_angles, 1);
% farfield_theta = farfield_theta + pi/2;
% farfield_phi = (pi/180)*[-180:1:179].';
% const_theta = 90; % used for plot
% % calculate farfields
% calculate_farfields_per_plane;
% %plotting the farfield data
% f = figure;
% pat1 = farfield_dataTheta(1,:).';
% pat2 = farfield_dataPhi(1,:).';
% 
% % if scale_type is db use these, otherwise comment these two lines
% pat1 = 10*log10(pat1);
% pat2 = 10*log10(pat2);
% 
% max_val = max(max([pat1 pat2]));
% max_val = step_size * ceil(max_val/step_size);
% 
% legend_str1 = ...
%     [plot_type '_{\theta}, f=' num2str(freq*1e-9) ' GHz'];
% legend_str2 = ...
%     [plot_type '_{\phi}, f=' num2str(freq*1e-9) ' GHz'];
% 
% polar_plot_constant_theta(farfield_phi,pat1,pat2,max_val, ...
%     step_size, Nrings,line_style1,line_style2,const_theta, ...
%     legend_str1,legend_str2,scale_type);

%% xz plane
% ===============================================
farfield_theta = zeros(number_of_angles, 1);
farfield_phi   = zeros(number_of_angles, 1);
farfield_theta = (pi/180)*[-179:1:0].';
const_phi = 0; % used for plot

% calculate farfields
calculate_farfields_per_plane;
% plotting the farfield data
 pat1 = farfield_dataTheta(1,:).';
 pat2 = farfield_dataPhi(1,:).';
% if scale_type is db use these, otherwise comment these two lines
 pat1 = 10*log10(pat1); 
 pat2 = 10*log10(pat2);
% 极化方向图的画法
%  max_val = max(max([pat1 pat2]));
%  max_val = step_size * ceil(max_val/step_size);
%  f = figure;
%  legend_str1 = ...
%  [plot_type '_{\theta}, f=' num2str(freq*1e-9) ' GHz'];
%  legend_str2 = ...
%  [plot_type '_{\phi}, f=' num2str(freq*1e-9) ' GHz'];
% 
%  polar_plot_constant_phi(farfield_theta,pat1,pat2,max_val, ...
%         step_size, Nrings,line_style1,line_style2,const_phi, ...
%         legend_str1,legend_str2,scale_type);

% 普通的plot

x = -179:0;
figure;plot(x,pat1);xlabel("bistatic angle/degree");ylabel("RCS/dB");title("HH polarzation");
figure;plot(x,pat2);xlabel("bistatic angle/degree");ylabel("RCS/dB");title("HV polarzation");
%% yz plane
% ===============================================
% farfield_theta = zeros(number_of_angles, 1);
% farfield_phi   = zeros(number_of_angles, 1);
% farfield_phi = farfield_phi + pi/2;
% farfield_theta = (pi/180)*[-180:1:179].';
% const_phi = 90; % used for plot
% 
% % calculate farfields
% calculate_farfields_per_plane;
% 
% % plotting the farfield data
%  f = figure;
%  pat1 = farfield_dataTheta(1,:).';
%  pat2 = farfield_dataPhi(1,:).';
% 
% % if scale_type is db use these, otherwise comment these two lines
%  pat1 = 10*log10(pat1); 
%  pat2 = 10*log10(pat2);
% 
%  max_val = max(max([pat1 pat2]));
%  max_val = step_size * ceil(max_val/step_size);
%  
%  legend_str1 = ...
%  [plot_type '_{\theta}, f=' num2str(freq*1e-9) ' GHz'];
%  legend_str2 = ...
%  [plot_type '_{\phi}, f=' num2str(freq*1e-9) ' GHz'];
% 
%  polar_plot_constant_phi(farfield_theta,pat1,pat2,max_val, ...
%         step_size, Nrings,line_style1,line_style2,const_phi, ...
%         legend_str1,legend_str2,scale_type);

