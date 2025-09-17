exp_jk_rpr = zeros(number_of_angles,1);
% dx_sinth_cosphi = zeros(number_of_angles,1);
% dy_sinth_sinphi = zeros(number_of_angles,1);
% dz_costh = zeros(number_of_angles,1);
% dy_dz_costh_sinphi = zeros(number_of_angles,1);
% dy_dz_sinth = zeros(number_of_angles,1);
% dy_dz_cosphi = zeros(number_of_angles,1);
% dx_dz_costh_cosphi = zeros(number_of_angles,1);
% dx_dz_sinth = zeros(number_of_angles,1);
% dx_dz_sinphi = zeros(number_of_angles,1);
% dx_dy_costh_cosphi = zeros(number_of_angles,1);
% dx_dy_costh_sinphi = zeros(number_of_angles,1);
% dx_dy_sinphi = zeros(number_of_angles,1);
% dx_dy_cosphi = zeros(number_of_angles,1);

farfield_dirTheta = zeros(1,number_of_angles);
farfield_dir = zeros(1,number_of_angles);
farfield_dirPhi = zeros(1,number_of_angles);

dx_sinth_cosphi = dx*sin(farfield_theta).*cos(farfield_phi);
dy_sinth_sinphi = dy*sin(farfield_theta).*sin(farfield_phi);
dz_costh = dz*cos(farfield_theta);
dy_dz_costh_sinphi = dy*dz*cos(farfield_theta).*sin(farfield_phi);
dy_dz_sinth = dy*dz*sin(farfield_theta);
dy_dz_cosphi = dy*dz*cos(farfield_phi);
dx_dz_costh_cosphi = dx*dz*cos(farfield_theta).*cos(farfield_phi);
dx_dz_sinth = dx*dz*sin(farfield_theta);
dx_dz_sinphi = dx*dz*sin(farfield_phi);
dx_dy_costh_cosphi = dx*dy*cos(farfield_theta).*cos(farfield_phi);
dx_dy_costh_sinphi = dx*dy*cos(farfield_theta).*sin(farfield_phi);
dx_dy_sinphi = dx*dy*sin(farfield_phi);
dx_dy_cosphi = dx*dy*cos(farfield_phi);
ci = 0.5*(ui+li);
cj = 0.5*(uj+lj);
ck = 0.5*(uk+lk);

% calculate directivity
k = 2*pi*freq*(mu_0*eps_0)^0.5;

Ntheta = zeros(number_of_angles,1);
Ltheta = zeros(number_of_angles,1);
Nphi = zeros(number_of_angles,1);
Lphi = zeros(number_of_angles,1);
rpr = zeros(number_of_angles,1);

for nj = lj:uj-1
    for nk =lk:uk-1
        % for +ax direction
        
        rpr = (ui - ci)*dx_sinth_cosphi ...
            + (nj-cj+0.5)*dy_sinth_sinphi ...
            + (nk-ck+0.5)*dz_costh;
        exp_jk_rpr = exp(j*k*rpr);
        Ntheta = Ntheta + (jyxp(1,nj-lj+1,nk-lk+1).*dy_dz_costh_sinphi ...
            - jzxp(1,nj-lj+1,nk-lk+1).*dy_dz_sinth).*exp_jk_rpr;
        Ltheta = Ltheta + (myxp(1,nj-lj+1,nk-lk+1).*dy_dz_costh_sinphi ...
            - mzxp(1,nj-lj+1,nk-lk+1).*dy_dz_sinth).*exp_jk_rpr;
        a1 = jyxp(1,nj-lj+1,nk-lk+1).*dy_dz_costh_sinphi;
        
        Nphi = Nphi + (jyxp(1,nj-lj+1,nk-lk+1).*dy_dz_cosphi).*exp_jk_rpr;
        Lphi = Lphi + (myxp(1,nj-lj+1,nk-lk+1).*dy_dz_cosphi).*exp_jk_rpr;
        
        % for -ax direction
        rpr = (li - ci)*dx_sinth_cosphi ...
            + (nj-cj+0.5)*dy_sinth_sinphi ...
            + (nk-ck+0.5)*dz_costh;
        exp_jk_rpr = exp(j*k*rpr);
        Ntheta = Ntheta + (jyxn(1,nj-lj+1,nk-lk+1).*dy_dz_costh_sinphi ...
            - jzxn(1,nj-lj+1,nk-lk+1).*dy_dz_sinth).*exp_jk_rpr;
        Ltheta = Ltheta + (myxn(1,nj-lj+1,nk-lk+1).*dy_dz_costh_sinphi ...
            - mzxn(1,nj-lj+1,nk-lk+1).*dy_dz_sinth).*exp_jk_rpr;
        Nphi = Nphi + (jyxn(1,nj-lj+1,nk-lk+1).*dy_dz_cosphi).*exp_jk_rpr;
        Lphi = Lphi + (myxn(1,nj-lj+1,nk-lk+1).*dy_dz_cosphi).*exp_jk_rpr;
    end
end
for ni =li:ui-1
    for nk =lk:uk-1
        % for +ay direction
        
        rpr = (ni - ci + 0.5)*dx_sinth_cosphi ...
            + (uj-cj)*dy_sinth_sinphi ...
            + (nk-ck+0.5)*dz_costh;
        exp_jk_rpr = exp(j*k*rpr);
        
        Ntheta = Ntheta + (jxyp(ni-li+1,1,nk-lk+1).*dx_dz_costh_cosphi ...
            - jzyp(ni-li+1,1,nk-lk+1).*dx_dz_sinth).*exp_jk_rpr;
        Ltheta = Ltheta + (mxyp(ni-li+1,1,nk-lk+1).*dx_dz_costh_cosphi ...
            - mzyp(ni-li+1,1,nk-lk+1).*dx_dz_sinth).*exp_jk_rpr;
        Nphi = Nphi + (-jxyp(ni-li+1,1,nk-lk+1).*dx_dz_sinphi).*exp_jk_rpr;
        Lphi = Lphi + (-mxyp(ni-li+1,1,nk-lk+1).*dx_dz_sinphi).*exp_jk_rpr;
        
        % for -ay direction
        rpr = (ni - ci + 0.5)*dx_sinth_cosphi ...
            + (lj-cj)*dy_sinth_sinphi ...
            + (nk-ck+0.5)*dz_costh;
        exp_jk_rpr = exp(j*k*rpr);
        
        Ntheta = Ntheta + (jxyn(ni-li+1,1,nk-lk+1).*dx_dz_costh_cosphi ...
            - jzyn(ni-li+1,1,nk-lk+1).*dx_dz_sinth).*exp_jk_rpr;
        Ltheta = Ltheta + (mxyn(ni-li+1,1,nk-lk+1).*dx_dz_costh_cosphi ...
            - mzyn(ni-li+1,1,nk-lk+1).*dx_dz_sinth).*exp_jk_rpr;
        Nphi = Nphi + (-jxyn(ni-li+1,1,nk-lk+1).*dx_dz_sinphi).*exp_jk_rpr;
        Lphi = Lphi + (-mxyn(ni-li+1,1,nk-lk+1).*dx_dz_sinphi).*exp_jk_rpr;
    end
end

for ni =li:ui-1
    for nj =lj:uj-1
        % for +az direction
        
        rpr = (ni-ci+0.5)*dx_sinth_cosphi ...
            + (nj - cj + 0.5)*dy_sinth_sinphi ...
            + (uk-ck)*dz_costh;
        exp_jk_rpr = exp(j*k*rpr);
        
        Ntheta = Ntheta + (jxzp(ni-li+1,nj-lj+1,1).*dx_dy_costh_cosphi ...
            + jyzp(ni-li+1,nj-lj+1,1).*dx_dy_costh_sinphi).*exp_jk_rpr;
        Ltheta = Ltheta + (mxzp(ni-li+1,nj-lj+1,1).*dx_dy_costh_cosphi ...
            + myzp(ni-li+1,nj-lj+1,1).*dx_dy_costh_sinphi).*exp_jk_rpr;
        Nphi = Nphi + (-jxzp(ni-li+1,nj-lj+1,1) ...
            .*dx_dy_sinphi+jyzp(ni-li+1,nj-lj+1,1).*dx_dy_cosphi).*exp_jk_rpr;
        Lphi = Lphi + (-mxzp(ni-li+1,nj-lj+1,1) ...
            .*dx_dy_sinphi+myzp(ni-li+1,nj-lj+1,1).*dx_dy_cosphi).*exp_jk_rpr;
        
        % for -az direction
        
        rpr = (ni-ci+0.5)*dx_sinth_cosphi ...
            + (nj - cj + 0.5)*dy_sinth_sinphi ...
            + (lk-ck)*dz_costh;
        exp_jk_rpr = exp(j*k*rpr);
        
        Ntheta = Ntheta + (jxzn(ni-li+1,nj-lj+1,1).*dx_dy_costh_cosphi ...
            + jyzn(ni-li+1,nj-lj+1,1).*dx_dy_costh_sinphi).*exp_jk_rpr;
        Ltheta = Ltheta + (mxzn(ni-li+1,nj-lj+1,1).*dx_dy_costh_cosphi ...
            + myzn(ni-li+1,nj-lj+1,1).*dx_dy_costh_sinphi).*exp_jk_rpr;
        Nphi = Nphi + (-jxzn(ni-li+1,nj-lj+1,1) ...
            .*dx_dy_sinphi+jyzn(ni-li+1,nj-lj+1,1).*dx_dy_cosphi).*exp_jk_rpr;
        Lphi = Lphi + (-mxzn(ni-li+1,nj-lj+1,1) ...
            .*dx_dy_sinphi+myzn(ni-li+1,nj-lj+1,1).*dx_dy_cosphi).*exp_jk_rpr;
    end
end

    % calculate radar cross-section
    farfield_dataTheta(1,:)  = ...
        (k^2./(8*pi*DEV.eta_0*radiated_power)) ...
        .* (abs(Lphi+DEV.eta_0*Ntheta).^2);
    farfield_dataPhi(1,:)    = ...
        (k^2./(8*pi*DEV.eta_0*radiated_power)) ...
        .* (abs(Ltheta-DEV.eta_0*Nphi).^2);
 
