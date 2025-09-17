uu = [];
for ii = 1:size(material_3d_space,1)
    for jj = 1:size(material_3d_space,2)
        for kk = 1:size(material_3d_space,3)
            if(material_3d_space(ii,jj,kk)~=1)
                uu = [uu;ii jj kk];
            end
        end
    end
end
figure;
x=uu(:,1); y=uu(:,2); z=uu(:,3); plot3(x,y,z,'o') ;

clear uu x y z