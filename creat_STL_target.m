cx = fdtd_domain.cell_center_coordinates_x;
cy = fdtd_domain.cell_center_coordinates_y;
cz = fdtd_domain.cell_center_coordinates_z;
n_ver = size(STL_target_x,2);
x_can = [];
for kk = 1:nz
    kk
    for jj = 1:ny
        y_c = cy(1,jj,1);
        z_c = cz(1,1,kk);
        for ind = 1:n_ver
            if any(STL_target_y(:,ind)>y_c) & any(STL_target_y(:,ind)<=y_c) & any(STL_target_z(:,ind)>z_c) & any(STL_target_z(:,ind)<=z_c)
                x_can = [x_can sum(STL_target_x(:,ind))/3];
            end
        end
        %plan A
        if ~isempty(x_can)
            max_x = max(x_can);
            min_x = min(x_can);
            max_x_idx = find_nearest(cx(:,jj,kk),max_x);
            min_x_idx = find_nearest(cx(:,jj,kk),min_x);
            for ii = min_x_idx:max_x_idx
                material_3d_space(ii,jj,kk) = STL_target.material_type;
            end
            x_can = [];
        end
    end
end

clear cx cy cz STL_target_x STL_target_y STL_target_z

function y = find_nearest(A,b)
[~,index]=sort(abs(A(:)-b));
y = index(1);
end

