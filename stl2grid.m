%%
filename = 'D:\Documents\Electromagnetic\model\空\Simplified Airplane\airplane_1.5G.txt';
C = load(filename);
fmax = 1.5e9;
c = 3e8;
dx = c/fmax/15;
dy=dx;dz=dx;
%%
% xyz方向上分别有多少个点
nx = (max(C(:,1)) - min(C(:,1)))/dx + 1;
ny = (max(C(:,2)) - min(C(:,2)))/dy + 1;
nz = (max(C(:,3)) - min(C(:,3)))/dz + 1;
nx = round(nx);
ny = round(ny);
nz = round(nz);

material_3d_space_temp = ones(nx,ny,nz);

% xyz方向上点的坐标（整数，从0开始，0,1,2,3,）
C_int(:,1) = (C(:,1) - min(C(:,1)))/dx;
C_int(:,2) = (C(:,2) - min(C(:,2)))/dy;
C_int(:,3) = (C(:,3) - min(C(:,3)))/dz;
C_int = round(C_int);

%每个点对应的material_3d_space元素的标号
C_I = C_int(:,3)*nx*ny + C_int(:,2)*nx + C_int(:,1) + 1;

% 将这些点赋值为金属
material_3d_space_temp(C_I) = 2;

% 补充空气层和CPML
material_3d_space = ones(nx+36,ny+36,nz+36);
material_3d_space(19:(19+nx-1),19:(19+ny-1),19:(19+nz-1)) = material_3d_space_temp;
[nx,ny,nz] = size(material_3d_space);

% 计算size
fdtd_domain.size_x = nx*dx;
fdtd_domain.size_y = ny*dy;
fdtd_domain.size_z = nz*dz;

fdtd_domain.min_x = -fdtd_domain.size_x/2;
fdtd_domain.min_y = -fdtd_domain.size_y/2;
fdtd_domain.min_z = 0;
fdtd_domain.max_x = fdtd_domain.size_x/2;
fdtd_domain.max_y = fdtd_domain.size_y/2;
fdtd_domain.max_z = fdtd_domain.size_z;

%%
% uu = [];
% for ii = 1:size(material_3d_space,1)
%     for jj = 1:size(material_3d_space,2)
%         for kk = 1:size(material_3d_space,3)
%             if(material_3d_space(ii,jj,kk)~=1)
%                 uu = [uu;ii jj kk];
%             end
%         end
%     end
% end
% figure;scatter3(uu(:,1),uu(:,2),uu(:,3),1);axis equal
% axis([1 nx 1 ny 1 nz]); 

%%
uu = [];
aaa = find(material_3d_space==2);
aaa = aaa-1;
uu(:,3) = floor(aaa/nx/ny);
aaa = aaa - uu(:,3)*nx*ny;
uu(:,2) = floor(aaa/nx);
uu(:,1) = aaa - uu(:,2)*nx;
uu = uu + 1;
figure;scatter3(uu(:,1),uu(:,2),uu(:,3),1);axis equal
axis([1 nx 1 ny 1 nz]); 

%%
filename_temp = filename(1:length(filename)-4);
filename_1 = strcat(filename_temp,'.txt');
filename_2 = strcat(filename_temp,'_ini.txt');

fid=fopen(filename_1,'w');
fprintf(fid,'%d ',material_3d_space);
fclose(fid);

fid=fopen(filename_2,'w');
fprintf(fid,'%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %d %d %d',...
    fdtd_domain.min_x,fdtd_domain.min_y,fdtd_domain.min_z,...
    fdtd_domain.max_x,fdtd_domain.max_y,fdtd_domain.max_z,...
    fdtd_domain.size_x,fdtd_domain.size_y,fdtd_domain.size_z,...
    nx,ny,nz);
fclose(fid);