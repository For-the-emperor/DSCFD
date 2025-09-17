figure;

for ind =1:size(material_types,2)
    matcol(ind,:) = material_types(ind).color;
end

m3d  = ones(nx+2, ny+2, nz+2);
tm3d = zeros(nx+2, ny+2, nz+2);
m3d(2:nx+1, 2:ny+1, 2:nz+1) = material_3d_space;
tm3d(2:nx+1, 2:ny+1, 2:nz+1) = m3d(2:nx+1, 2:ny+1, 2:nz+1) ...
    - m3d(1:nx, 2:ny+1, 2:nz+1);
I1 = find(tm3d>0);
tm3d(2:nx+1, 2:ny+1, 2:nz+1) = m3d(2:nx+1, 2:ny+1, 2:nz+1) ...
    - m3d(3:nx+2, 2:ny+1, 2:nz+1);
I2 = find(tm3d>2);
tm3d(2:nx+1, 2:ny+1, 2:nz+1) = m3d(2:nx+1, 2:ny+1, 2:nz+1) ...
    - m3d(2:nx+1, 1:ny, 2:nz+1);
I3 = find(tm3d>2);
tm3d(2:nx+1, 2:ny+1, 2:nz+1) = m3d(2:nx+1, 2:ny+1, 2:nz+1) ...
    - m3d(2:nx+1, 3:ny+2, 2:nz+1);
I4 = find(tm3d>2);
tm3d(2:nx+1, 2:ny+1, 2:nz+1) = m3d(2:nx+1, 2:ny+1, 2:nz+1) ...
    - m3d(2:nx+1, 2:ny+1, 1:nz);
I5 = find(tm3d>2);
tm3d(2:nx+1, 2:ny+1, 2:nz+1) = m3d(2:nx+1, 2:ny+1, 2:nz+1) ...
    - m3d(2:nx+1, 2:ny+1, 3:nz+2);
I6 = find(tm3d>2);
clear tm3d;

xx = zeros(nx+2,ny+2,nz+2);
yy = zeros(nx+2,ny+2,nz+2);
zz = zeros(nx+2,ny+2,nz+2);

for ind = -1:nx
    xx(ind+2,:,:) = ind * dx;
end
xx = xx + fdtd_domain.min_x;
for ind = -1:ny
    yy(:,ind+2,:) = ind * dy;
end
yy = yy + fdtd_domain.min_y;
for ind = -1:nz
    zz(:,:,ind+2) = ind * dz;
end
zz = zz + fdtd_domain.min_z;

len = size(I1,1);
vx = [xx(I1) xx(I1) xx(I1) xx(I1)];
vy = [yy(I1) yy(I1)+dy  yy(I1)+dy yy(I1)];
vz = [zz(I1) zz(I1) zz(I1)+dz zz(I1)+dz];
vx  = reshape(vx.',4*len,1);    
vy  = reshape(vy.',4*len,1);    
vz  = reshape(vz.',4*len,1);    
v1 = [vx vy vz];
f1 = [(1:4:4*(len-1)+1).' (2:4:4*(len-1)+2).' (3:4:4*(len-1)+3).' (4:4:4*(len-1)+4).'];
RGB1 = matcol(m3d(I1),:);

len = size(I2,1);
vx = [xx(I2) xx(I2) xx(I2) xx(I2)]+dx;
vy = [yy(I2) yy(I2)+dy  yy(I2)+dy yy(I2)];
vz = [zz(I2) zz(I2) zz(I2)+dz zz(I2)+dz];
vx  = reshape(vx.',4*len,1);    
vy  = reshape(vy.',4*len,1);    
vz  = reshape(vz.',4*len,1);    
v2 = [vx vy vz];
f2 = [(1:4:4*(len-1)+1).' (2:4:4*(len-1)+2).' (3:4:4*(len-1)+3).' (4:4:4*(len-1)+4).'];
RGB2 = matcol(m3d(I2),:);

len = size(I3,1);
vx = [xx(I3) xx(I3)+dx xx(I3)+dx xx(I3)];
vy = [yy(I3) yy(I3) yy(I3) yy(I3)];
vz = [zz(I3) zz(I3) zz(I3)+dz zz(I3)+dz];
vx  = reshape(vx.',4*len,1);    
vy  = reshape(vy.',4*len,1);    
vz  = reshape(vz.',4*len,1);    
v3 = [vx vy vz];
f3 = [(1:4:4*(len-1)+1).' (2:4:4*(len-1)+2).' (3:4:4*(len-1)+3).' (4:4:4*(len-1)+4).'];
RGB3 = matcol(m3d(I3),:);

len = size(I4,1);
vx = [xx(I4) xx(I4)+dx xx(I4)+dx xx(I4)];
vy = [yy(I4) yy(I4) yy(I4) yy(I4)] + dy;
vz = [zz(I4) zz(I4) zz(I4)+dz zz(I4)+dz];
vx  = reshape(vx.',4*len,1);    
vy  = reshape(vy.',4*len,1);    
vz  = reshape(vz.',4*len,1);    
v4 = [vx vy vz];
f4 = [(1:4:4*(len-1)+1).' (2:4:4*(len-1)+2).' (3:4:4*(len-1)+3).' (4:4:4*(len-1)+4).'];
RGB4 = matcol(m3d(I4),:);

len = size(I5,1);
vx = [xx(I5) xx(I5) xx(I5)+dx xx(I5)+dx];
vy = [yy(I5) yy(I5)+dy  yy(I5)+dy yy(I5)];
vz = [zz(I5) zz(I5) zz(I5) zz(I5)];
vx  = reshape(vx.',4*len,1);    
vy  = reshape(vy.',4*len,1);    
vz  = reshape(vz.',4*len,1);    
v5 = [vx vy vz];
f5 = [(1:4:4*(len-1)+1).' (2:4:4*(len-1)+2).' (3:4:4*(len-1)+3).' (4:4:4*(len-1)+4).'];
RGB5 = matcol(m3d(I5),:);

len = size(I6,1);
vx = [xx(I6) xx(I6) xx(I6)+dx xx(I6)+dx];
vy = [yy(I6) yy(I6)+dy  yy(I6)+dy yy(I6)];
vz = [zz(I6) zz(I6) zz(I6) zz(I6)] + dz;
vx  = reshape(vx.',4*len,1);    
vy  = reshape(vy.',4*len,1);    
vz  = reshape(vz.',4*len,1);    
v6 = [vx vy vz];
f6 = [(1:4:4*(len-1)+1).' (2:4:4*(len-1)+2).' (3:4:4*(len-1)+3).' (4:4:4*(len-1)+4).'];
RGB6 = matcol(m3d(I6),:);

vertices = [v1;v2;v3;v4;v5;v6];

faces = f1; 
faces = [faces; f2+ max(max(faces))];
faces = [faces; f3+ max(max(faces))];
faces = [faces; f4+ max(max(faces))];
faces = [faces; f5+ max(max(faces))];
faces = [faces; f6+ max(max(faces))];
RGB  = [RGB1;RGB2;RGB3;RGB4;RGB5;RGB6];

if ~isempty(vertices)
    patch('Vertices',vertices,'Faces',faces,'facevertexcdata',RGB,'FaceColor','flat');
end

[i,j,k]=ind2sub(size(material_3d_space),find(material_3d_space~=1));
figure;
scatter3(i,j,k);