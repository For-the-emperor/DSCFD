%% 用跟FDTD类似地方法写
%% 定义简单的目标结构
disp("defining the targets");

bricks = [];
spheres = [];

% 定义立方体，采用六个边界的坐标
% bricks(1).min_x = -0.25;
% bricks(1).min_y = -0.25;
% bricks(1).min_z = 0;
% bricks(1).max_x = 0.25;
% bricks(1).max_y = 0.25;
% bricks(1).max_z = 0.5;
% bricks(1).material_type = 2 ;

%定义球体，采用球的半径和圆心坐标
spheres(1).radius = 0.9;
spheres(1).center_x = 0;
spheres(1).center_y = 0;
spheres(1).center_z = 0;
spheres(1).material_type = 2 ;

% 计算计算域的尺寸
disp("calculating domain size")
number_of_bricks = size(bricks,2);
numbber_of_spheres = size(spheres,2);

number_of_objects = 1;
for i=1:numbber_of_spheres
    min_x(number_of_objects) = spheres(i).center_x - spheres(i).radius; 
    min_y(number_of_objects) = spheres(i).center_y - spheres(i).radius; 
    min_z(number_of_objects) = spheres(i).center_z - spheres(i).radius; 
    max_x(number_of_objects) = spheres(i).center_x + spheres(i).radius; 
    max_y(number_of_objects) = spheres(i).center_y + spheres(i).radius; 
    max_z(number_of_objects) = spheres(i).center_z + spheres(i).radius; 
    number_of_objects = number_of_objects + 1;
end

for i=1:number_of_bricks
    min_x(number_of_objects) = bricks(i).min_x; 
    min_y(number_of_objects) = bricks(i).min_y; 
    min_z(number_of_objects) = bricks(i).min_z; 
    max_x(number_of_objects) = bricks(i).max_x; 
    max_y(number_of_objects) = bricks(i).max_y; 
    max_z(number_of_objects) = bricks(i).max_z; 
    number_of_objects = number_of_objects + 1;
end
% 引入空气层和CPML边界后，计算区域尺寸
min_x = min(min_x) - dx*DEV.air_buffer - dx*DEV.NPML(1);
min_y = min(min_y) - dy*DEV.air_buffer - dy*DEV.NPML(3);
min_z = min(min_z) - dz*DEV.air_buffer - dz*DEV.NPML(5);
max_x = max(max_x) + dx*DEV.air_buffer + dx*DEV.NPML(2);
max_y = max(max_y) + dy*DEV.air_buffer + dy*DEV.NPML(4);
max_z = max(max_z) + dz*DEV.air_buffer + dz*DEV.NPML(6);

Sx = max_x-min_x; Sy = max_y - min_y; Sz = max_z - min_z;
% 调整计算区域，使其为网格尺寸的整数倍
Nx = ceil(Sx/dx);
Sx = Nx*dx;
max_x = min_x + Sx;
Ny = ceil(Sy/dy);
Sy = Ny*dy;
max_y = min_y + Sy;
Nz = ceil(Sz/dz);
Sz = Nz*dz;
max_z = min_z + Sz;
% 2x GRID
Nx2 = 2*Nx; dx2 = dx/2;
Ny2 = 2*Ny; dy2 = dy/2;
Nz2 = 2*Nz; dz2 = dz/2;

% 网格中心坐标
xa = (min_x+dx/2:dx:max_x);
ya = (min_y+dy/2:dy:max_y);
za = (min_z+dz/2:dz:max_z);
% 生成存储三个坐标的三维矩阵
[Y,X,Z] = meshgrid(ya,xa,za);

% 双倍网格同理
xa2 = (min_x+dx2/2:dx2:max_x);
ya2 = (min_y+dy2/2:dy2:max_y);
za2 = (min_z+dz2/2:dz2:max_z);
[Y2,X2,Z2] = meshgrid(ya2,xa2,za2);

%% 电磁参数赋值
disp("initializing target");
% 初始化
DEV.ER2 = ones(Nx2,Ny2,Nz2)*material_types(1).eps;
DEV.UR2 = ones(Nx2,Ny2,Nz2)*material_types(1).mu;
% 建立目标
% 球
for ind = 1:numbber_of_spheres
    distance = sqrt((spheres(ind).center_x - X2).^2 + (spheres(ind).center_y - Y2).^2 + (spheres(ind).center_z - Z2).^2);
    I = find(distance<=spheres(ind).radius);
    DEV.ER2(I) = material_types(spheres(ind).material_type).eps;
    DEV.UR2(I) = material_types(spheres(ind).material_type).mu;
end

% 立方体
for ind = 1:number_of_bricks
    blx = round((bricks(ind).min_x - min_x)/dx2)+1;
    bly = round((bricks(ind).min_y - min_y)/dy2)+1;
    blz = round((bricks(ind).min_z - min_z)/dz2)+1;

    bux = round((bricks(ind).max_x - min_x)/dx2)+1;
    buy = round((bricks(ind).max_y - min_y)/dy2)+1;
    buz = round((bricks(ind).max_z - min_z)/dz2)+1;

    DEV.ER2(blx:bux-1,bly:buy-1,blz:buz-1) = material_types(bricks(ind).material_type).eps;
    DEV.UR2(blx:bux-1,bly:buy-1,blz:buz-1) = material_types(bricks(ind).material_type).mu;

end