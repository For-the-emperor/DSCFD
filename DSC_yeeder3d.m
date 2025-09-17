function [DEX,DEY,DEZ,DHX,DHY,DHZ] = DSC_yeeder3d(NS,RES,BC,kinc,bandwidth)
% YEEDER3D      Derivative Matrices on a 3D Yee Grid
%
% [DEX,DEY,DEZ,DHX,DHY,DHZ] = yeeder3d(NS,RES,BC,kinc);
%
% INPUT ARGUMENTS
% =================
% NS    [Nx Ny Nz] Grid Size
% RES   [dx dy dz] Grid Resolution
% BC    [xbc ybc zbc] Boundary Conditions
%         0: Dirichlet boundary conditions
%         1: Periodic boundary conditions
% kinc  [kx ky kz] Incident Wave Vector
%       This argument is only needed for PBCs.
% M     Bandwidth of DSC scheme (2,3 or 4)
%
% Note: For normalized grids use k0*RES and kinc/k0
%
% OUTPUT ARGUMENTS
% =================
% DEX   Derivative Matrix wrt x for Electric Fields
% DEY   Derivative Matrix wrt y for Electric Fields
% DEZ   Derivative Matrix wrt z for Electric Fields
% DHX   Derivative Matrix wrt x for Magnetic Fields
% DHY   Derivative Matrix wrt y for Magnetic Fields
% DHZ   Derivative Matrix wrt z for Magnetic Fields

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HANDLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EXTRACT GRID PARAMETERS
Nx = NS(1);     dx = RES(1);
Ny = NS(2);     dy = RES(2);
Nz = NS(3);     dz = RES(3);

% DEFAULT KINC
if ~exist('kinc')
    kinc = [0 0 0];
end

% DETERMINE MATRIX SIZE
M = Nx*Ny*Nz;

% ZERO MATRIX
Z = sparse(M,M);

%define DSC coefficients
if bandwidth == 1
    coeffs = [-1, 1];
    offsets = [0, 1];
    min_N = 3;
    
elseif bandwidth == 2
    coeffs = [1/24, -9/8, 9/8, -1/24];
    offsets = [-1, 0, 1, 2];
    min_N = 5;
elseif bandwidth == 3
    coeffs = [-3/640, 25/384, -75/64, 75/64, -25/384, 3/640];
    offsets = [-2, -1, 0, 1, 2, 3];
    min_N = 7;
elseif bandwidth == 4
    coeffs = [5/7168, -49/5120, 245/3072, -1225/1024, 1225/1024, -245/3072, 49/5120, -5/7168]; % x-4 to x+4
    offsets = [-3, -2,-1, 0, 1, 2, 3, 4];
    min_N = 9;
else 
    error('Unsupported M value. Use M=2, 3, or 4.');
end
num_of_diags = length(offsets);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 为了方便施加边界条件，给每一行提供xyz的索引
idx3 = (1:M).';               % 列向量
ix   = mod(idx3-1, Nx) + 1;                       % x 从 1 到 Nx 循环
iy   = mod(floor((idx3-1)/Nx), Ny) + 1;           % y 从 1 到 Ny 循环
iz   = floor((idx3-1)/(Nx*Ny)) + 1;               % z 从 1 到 Nz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DEX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Nx < min_N
    error(['Grid size Nx too small for M=', num2str(M)]);
end
diags = zeros(M,num_of_diags);
for i = 1:num_of_diags
    diags(:,i) = coeffs(i);
end

for i = 1:num_of_diags
    o = offsets(i);
    if o == 0
        continue;   % 主对角线永远合法
    end

    k = idx3;       % 源点索引向量
    j = k - o;      % offset后，diag所在的行号
    %包含两部分，1：超出计算域，2：主对角线与offset对角线不在同一个Y-Z平面上
%     invalid = (j < 1) | (j > M) | floor((j-1)/Nx) ~= floor((k-1)/Nx); 
    invalid = floor((j-1)/Nx) ~= floor((k-1)/Nx);
    diags(invalid, i) = 0;
end

DEX = spdiags(diags/dx, offsets, Z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DEY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Ny < min_N
    error(['Grid size Ny too small for M=', num2str(M)]);
end

% 按 y 方向偏移
diags = zeros(M,num_of_diags);
for i = 1:num_of_diags
    diags(:,i) = coeffs(i);
end
for i = 1:num_of_diags
    o = offsets(i);
    if o == 0
        continue;   % 主对角线永远合法
    end
    k = idx3;                 % 源点索引
    j = k - o * Nx;           % y 方向偏移 Nx 格
    % 越界 或 跨越 x-z 平面
%     invalid = (j < 1) | (j > M) | (floor((j-1)/(Nx*Ny)) ~= floor((k-1)/(Nx*Ny)));
    invalid = (floor((j-1)/(Nx*Ny)) ~= floor((k-1)/(Nx*Ny)));
    diags(invalid, i) = 0;
end
DEY = spdiags(diags / dy, offsets * Nx, Z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DEZ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Nz < min_N
    error(['Grid size Nz too small for M=', num2str(M)]);
end
% 按 z 方向偏移
diags = zeros(M,num_of_diags);
for i = 1:num_of_diags
    diags(:,i) = coeffs(i);
end
for i = 1:num_of_diags
    o = offsets(i);
    if o == 0
        continue;   % 主对角线永远合法
    end
    k = idx3;                   % 源点索引
    j = k - o * (Nx * Ny);      % z 方向偏移 Nx*Ny 格
    % 仅需检查越界
    invalid = (j < 1) | (j > M);
    diags(invalid, i) = 0;
end
DEZ = spdiags(diags / dz, offsets * Nx * Ny, Z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DHX, DHY AND DHZ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DHX = -DEX';
DHY = -DEY';
DHZ = -DEZ';