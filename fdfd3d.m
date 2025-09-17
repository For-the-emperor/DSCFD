function DAT = fdfd3d(DEV,SRC)
% FDFD3D 三维FDFD仿真程序
%
% DAT = fdfd3d(DEV,SRC);
%
% 输入参数
% ================
%
% .ER2 or
% .ER2xx .ER2xy .ER2xz
% .ER2yx .ER2yy .ER2yz relative perm. on 2x grid
% .ER2zx .ER2zy .ER2zz 两倍网格下的相对介电常数
%
% .UR2 or
% .UR2xx .UR2xy .UR2xz
% .UR2yx .UR2yy .UR2yz relative perm. on 2x grid
% .UR2zx .UR2zy .UR2zz 两倍网格下的磁导率
% 由于给Maxwell方程组做了归一化，完全不需要eps0和mu0
%
% .RES [dx dy] Yee网格尺寸
% .NPML [NZLO NZHI] CPML层数
%
% SRC 入射波参数
%
% .lam0 波长
% .theta 入射角theta，与Z轴正方向夹角
% .phi 入射角phi，与x轴正方向夹角
% .pte TE极化的复振幅
% .ptm TM极化的复振幅
%
% OUTPUT ARGUMENTS
% ================
% 输出
% DAT结构体
% .Ex .Ey .Ez 电场分量
% .Hx .Hy .Hz 磁场分量

% DEFINE SOLVER PARAMETERS
tol = 1e-5;
maxit = 15000;

% ANONYMOUS FUNCTIONS
diagonalize = @(x) diag(sparse(x(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HANDLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bandwidth = DEV.bandwidth;
% GET SIZE OF GRID
if isfield(DEV,'ER2')
    [Nx2,Ny2,Nz2] = size(DEV.ER2);
else
    [Nx2,Ny2,Nz2] = size(DEV.ER2xx);
end

% CALCULATE GRID
Nx = Nx2/2;
Ny = Ny2/2;
Nz = Nz2/2;
dx = DEV.RES(1);
dy = DEV.RES(2);
dz = DEV.RES(3);
Sx = Nx*dx;
Sy = Ny*dy;
Sz = Nz*dz;

% GRID AXES
xa2 = [1:Nx2]*dx/2; xa2 = xa2 - mean(xa2);
ya2 = [1:Ny2]*dy/2; ya2 = ya2 - mean(ya2);
za2 = [0.5:Nz2-0.5]*dz/2;
[Y2,X2,Z2] = meshgrid(ya2,xa2,za2);

% WAVE NUMBER
k0 = 2*pi/SRC.lam0;

% SPECIAL MATRICES
M = Nx*Ny*Nz;
ZZ = sparse(M,M);
I = speye(M,M);

% PERMITTIVITY TENSOR ELEMENTS
if isfield(DEV,'ER2')
    ERxx = diagonalize(DEV.ER2(2:2:Nx2,1:2:Ny2,1:2:Nz2));
    ERxy = ZZ;
    ERxz = ZZ;
    ERyx = ZZ;
    ERyy = diagonalize(DEV.ER2(1:2:Nx2,2:2:Ny2,1:2:Nz2));
    ERyz = ZZ;
    ERzx = ZZ;
    ERzy = ZZ;
    ERzz = diagonalize(DEV.ER2(1:2:Nx2,1:2:Ny2,2:2:Nz2));
else
    ERxx = diagonalize(DEV.ER2xx(2:2:Nx2,1:2:Ny2,1:2:Nz2));
    ERxy = diagonalize(DEV.ER2xy(1:2:Nx2,2:2:Ny2,1:2:Nz2));
    ERxz = diagonalize(DEV.ER2xz(1:2:Nx2,1:2:Ny2,2:2:Nz2));
    ERyx = diagonalize(DEV.ER2yx(2:2:Nx2,1:2:Ny2,1:2:Nz2));
    ERyy = diagonalize(DEV.ER2yy(1:2:Nx2,2:2:Ny2,1:2:Nz2));
    ERyz = diagonalize(DEV.ER2yz(1:2:Nx2,1:2:Ny2,2:2:Nz2));
    ERzx = diagonalize(DEV.ER2zx(2:2:Nx2,1:2:Ny2,1:2:Nz2));
    ERzy = diagonalize(DEV.ER2zy(1:2:Nx2,2:2:Ny2,1:2:Nz2));
    ERzz = diagonalize(DEV.ER2zz(1:2:Nx2,1:2:Ny2,2:2:Nz2));
end

% PERMEABILITY TENSOR ELEMENTS
if isfield(DEV,'UR2')
    URxx = diagonalize(DEV.UR2(1:2:Nx2,2:2:Ny2,2:2:Nz2));
    URxy = ZZ;
    URxz = ZZ;
    URyx = ZZ;
    URyy = diagonalize(DEV.UR2(2:2:Nx2,1:2:Ny2,2:2:Nz2));
    URyz = ZZ;
    URzx = ZZ;
    URzy = ZZ;
    URzz = diagonalize(DEV.UR2(2:2:Nx2,2:2:Ny2,1:2:Nz2));
else
    URxx = diagonalize(DEV.UR2xx(1:2:Nx2,2:2:Ny2,2:2:Nz2));
    URxy = diagonalize(DEV.UR2xy(2:2:Nx2,1:2:Ny2,2:2:Nz2));
    URxz = diagonalize(DEV.UR2xz(2:2:Nx2,2:2:Ny2,1:2:Nz2));
    URyx = diagonalize(DEV.UR2yx(1:2:Nx2,2:2:Ny2,2:2:Nz2));
    URyy = diagonalize(DEV.UR2yy(2:2:Nx2,1:2:Ny2,2:2:Nz2));
    URyz = diagonalize(DEV.UR2yz(2:2:Nx2,2:2:Ny2,1:2:Nz2));
    URzx = diagonalize(DEV.UR2zx(1:2:Nx2,2:2:Ny2,2:2:Nz2));
    URzy = diagonalize(DEV.UR2zy(2:2:Nx2,1:2:Ny2,2:2:Nz2));
    URzz = diagonalize(DEV.UR2zz(2:2:Nx2,2:2:Ny2,1:2:Nz2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM FDFD ANALYS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EXTRACT MATERIAL PROPERTIES IN EXTERNAL REGIONS
ur_ref = full(URxx(1,1));
ur_trn = full(URxx(M,M));
er_ref = full(ERxx(1,1));
er_trn = full(ERxx(M,M));
n_ref = sqrt(ur_ref*er_ref);
n_trn = sqrt(ur_trn*er_trn);

% CALCULATE PML PARAMETERS
[sx2,sy2,sz2] = calcpml3d([Nx2 Ny2 Nz2],2*DEV.NPML);

sx_ey = diagonalize(1./sx2(1:2:Nx2,2:2:Ny2,1:2:Nz2));
sx_ez = diagonalize(1./sx2(1:2:Nx2,1:2:Ny2,2:2:Nz2));
sy_ex = diagonalize(1./sy2(2:2:Nx2,1:2:Ny2,1:2:Nz2));
sy_ez = diagonalize(1./sy2(1:2:Nx2,1:2:Ny2,2:2:Nz2));
sz_ex = diagonalize(1./sz2(2:2:Nx2,1:2:Ny2,1:2:Nz2));
sz_ey = diagonalize(1./sz2(1:2:Nx2,2:2:Ny2,1:2:Nz2));

sx_hy = diagonalize(1./sx2(2:2:Nx2,1:2:Ny2,2:2:Nz2));
sx_hz = diagonalize(1./sx2(2:2:Nx2,2:2:Ny2,1:2:Nz2));
sy_hx = diagonalize(1./sy2(1:2:Nx2,2:2:Ny2,2:2:Nz2));
sy_hz = diagonalize(1./sy2(2:2:Nx2,2:2:Ny2,1:2:Nz2));
sz_hx = diagonalize(1./sz2(1:2:Nx2,2:2:Ny2,2:2:Nz2));
sz_hy = diagonalize(1./sz2(2:2:Nx2,1:2:Ny2,2:2:Nz2));

% CALCULATE WAVE VECTORS
kinc = k0*n_ref*[ sin(SRC.theta)*cos(SRC.phi) ; ...
    sin(SRC.theta)*sin(SRC.phi) ; ...
    cos(SRC.theta) ];
% m = [-floor(Nx/2):+floor((Nx-1)/2)].';
% n = [-floor(Ny/2):+floor((Ny-1)/2)].';
% kx = kinc(1) - m*2*pi/Sx;
% ky = kinc(2) - n*2*pi/Sy;
% kz_ref = sqrt((k0*n_ref)^2 - kx.^2 - ky.^2);
% kz_trn = sqrt((k0*n_trn)^2 - kx.^2 - ky.^2);

% BUILD DERIVATIVE MATRICES
NS = [Nx Ny Nz];
RES = [dx dy dz];
BC = [0 0 0];
[DEX,DEY,DEZ,DHX,DHY,DHZ] = DSC_yeeder3d(NS,k0*RES,BC,kinc/k0, bandwidth);

% CALCULATE INTERPOLATION MATRICES (Assumes SRC.theta = 0)
RX = (0.5*k0*dx)*abs(DEX);
RY = (0.5*k0*dy)*abs(DEY);
RZ = (0.5*k0*dz)*abs(DEZ);

% FORM MATERIALS TENSORS
ER = [ ERxx RX*RY'*ERxy RX*RZ'*ERxz ; ...
    RY*RX'*ERyx ERyy RY*RZ'*ERyz ; ...
    RZ*RX'*ERzx RZ*RY'*ERzy ERzz ];
UR = [ URxx RX'*RY*URxy RX'*RZ*URxz ; ...
    RY'*RX*URyx URyy RY'*RZ*URyz ; ...
    RZ'*RX*URzx RZ'*RY*URzy URzz ];

% BUILD THE WAVE MATRIX A
CE = [ ZZ , -sz_hx*DEZ , sy_hx*DEY ; ...
    sz_hy*DEZ , ZZ , -sx_hy*DEX ; ...
    -sy_hz*DEY , sx_hz*DEX , ZZ ];
CH = [ ZZ , -sz_ex*DHZ , sy_ex*DHY ; ...
    sz_ey*DHZ , ZZ , -sx_ey*DHX ; ...
    -sy_ez*DHY , sx_ez*DHX , ZZ ];
[i,j,v]=find(UR);
UR_inv = sparse(i,j,1./v,size(UR,1),size(UR,2));
A_e = CH*UR_inv*CE-ER;
% A_e = CH/UR*CE - ER;
% 清理变量释放内存
clear CH ER DEX DEY DEZ DHX DHY DHZ sx_ez sx_ey sx_hz sx_hy sy_ez sy_ex sy_hz sy_hx sz_ey sz_ex sz_hy sz_hx;
% CALCULATE POLARIZATION VECTOR P
az = [0;0;1];
if abs(SRC.theta)<1e-6
    ate = [0;1;0];
else
    ate = cross(az,kinc);
    ate = ate/norm(ate);
end
atm = cross(kinc,ate);
atm = atm/norm(atm);
P = SRC.pte*ate + SRC.ptm*atm;

% CALCULATE SOURCE FIELD fsrc
fphase = exp(-1i*(kinc(1)*X2 + kinc(2)*Y2 + kinc(3)*Z2));
fx = P(1)*fphase(2:2:Nx2,1:2:Ny2,1:2:Nz2);
fy = P(2)*fphase(1:2:Nx2,2:2:Ny2,1:2:Nz2);
fz = P(3)*fphase(1:2:Nx2,1:2:Ny2,2:2:Nz2);
fsrc = [ fx(:) ; fy(:) ; fz(:) ];

% CALCULATE SCATERED FIELD MASKING MATRIX
% 散射场：1；总场：0
n = DEV.TFSF;
Q = ones(Nx,Ny,Nz);
Q(n+1:Nx-n,n+1:Ny-n,n+1:Nz-n) = 0;
Q = diagonalize(Q);
Q = [ Q ZZ ZZ ; ZZ Q ZZ ; ZZ ZZ Q ];

% CALCULATE SOURCE VECTOR B
b = (Q*A_e - A_e*Q)*fsrc;

%% SOLVE ELECTRIC FIELD
f = zeros(3*M,1);

A_e = gpuArray(A_e);
b = gpuArray(b);
f = bicg(A_e,b,tol,maxit,[],[],f);
f = gather(f);


% EXTRACT FIELD COMPONENTS
DAT.Ex = reshape(f( 1: M),Nx,Ny,Nz);
DAT.Ey = reshape(f( M+1:2*M),Nx,Ny,Nz);
DAT.Ez = reshape(f(2*M+1:3*M),Nx,Ny,Nz);
clear A_e b;
%% SOLVE MAGNETIC FIELD
h = UR\CE*f;

DAT.Hx = reshape(h( 1: M),Nx,Ny,Nz);
DAT.Hy = reshape(h( M+1:2*M),Nx,Ny,Nz);
DAT.Hz = reshape(h(2*M+1:3*M),Nx,Ny,Nz);
DAT.Hx = 1i*DAT.Hx/DEV.eta_0;
DAT.Hy = 1i*DAT.Hy/DEV.eta_0;
DAT.Hz = 1i*DAT.Hz/DEV.eta_0;
