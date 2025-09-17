jxyp = zeros(ui-li,1,uk-lk);
jxzp = zeros(ui-li,uj-lj,1);
jyxp = zeros(1,uj-lj,uk-lk);
jyzp = zeros(ui-li,uj-lj,1);
jzxp = zeros(1,uj-lj,uk-lk);
jzyp = zeros(ui-li,1,uk-lk);
jxyn = zeros(ui-li,1,uk-lk);
jxzn = zeros(ui-li,uj-lj,1);
jyxn = zeros(1,uj-lj,uk-lk);
jyzn = zeros(ui-li,uj-lj,1);
jzxn = zeros(1,uj-lj,uk-lk);
jzyn = zeros(ui-li,1,uk-lk);

mxyp = zeros(ui-li,1,uk-lk);
mxzp = zeros(ui-li,uj-lj,1);
myxp = zeros(1,uj-lj,uk-lk);
myzp = zeros(ui-li,uj-lj,1);
mzxp = zeros(1,uj-lj,uk-lk);
mzyp = zeros(ui-li,1,uk-lk);
mxyn = zeros(ui-li,1,uk-lk);
mxzn = zeros(ui-li,uj-lj,1);
myxn = zeros(1,uj-lj,uk-lk);
myzn = zeros(ui-li,uj-lj,1);
mzxn = zeros(1,uj-lj,uk-lk);
mzyn = zeros(ui-li,1,uk-lk);

% 进行外推边界表面电磁流的计算
myxp(1,:,:) = 0.5 * (DAT.Ez(ui,lj:uj-1,lk:uk-1)+DAT.Ez(ui,lj+1:uj,lk:uk-1));
mzxp(1,:,:) = -0.5* (DAT.Ey(ui,lj:uj-1,lk:uk-1)+DAT.Ey(ui,lj:uj-1,lk+1:uk));
mxyp(:,1,:) = -0.5* (DAT.Ez(li:ui-1,uj,lk:uk-1)+DAT.Ez(li+1:ui,uj,lk:uk-1));
mzyp(:,1,:) = 0.5 * (DAT.Ex(li:ui-1,uj,lk:uk-1)+DAT.Ex(li:ui-1,uj,lk+1:uk));
mxzp(:,:,1) = 0.5 * (DAT.Ey(li:ui-1,lj:uj-1,uk)+DAT.Ey(li+1:ui,lj:uj-1,uk));
myzp(:,:,1) = -0.5* (DAT.Ex(li:ui-1,lj:uj-1,uk)+DAT.Ex(li:ui-1,lj+1:uj,uk));

jyxp(1,:,:) = -0.25*(DAT.Hz(ui,lj:uj-1,lk:uk-1)+DAT.Hz(ui,lj:uj-1,lk+1:uk)+DAT.Hz(ui-1,lj:uj-1,lk:uk-1)+DAT.Hz(ui-1,lj:uj-1,lk+1:uk));
jzxp(1,:,:) = 0.25 *(DAT.Hy(ui,lj:uj-1,lk:uk-1)+DAT.Hy(ui,lj+1:uj,lk:uk-1)+DAT.Hz(ui-1,lj:uj-1,lk:uk-1)+DAT.Hz(ui-1,lj+1:uj,lk:uk-1));
jzyp(:,1,:) = -0.25*(DAT.Hx(li:ui-1,uj,lk:uk-1)+DAT.Hx(li+1:ui,uj,lk:uk-1)+DAT.Hx(li:ui-1,uj-1,lk:uk-1)+DAT.Hx(li+1:ui,uj-1,lk:uk-1));
jxyp(:,1,:) = 0.25 *(DAT.Hz(li:ui-1,uj,lk:uk-1)+DAT.Hz(li:ui-1,uj,lk+1:uk)+DAT.Hz(li:ui-1,uj-1,lk:uk-1)+DAT.Hz(li:ui-1,uj-1,lk+1:uk));
jyzp(:,:,1) = 0.25 *(DAT.Hx(li:ui-1,lj:uj-1,uk)+DAT.Hx(li+1:ui,lj:uj-1,uk)+DAT.Hx(li:ui-1,lj:uj-1,uk-1)+DAT.Hx(li+1:ui,lj:uj-1,uk-1));
jxzp(:,:,1) = -0.25*(DAT.Hy(li:ui-1,lj:uj-1,uk)+DAT.Hy(li:ui-1,lj+1:uj,uk)+DAT.Hy(li:ui-1,lj:uj-1,uk-1)+DAT.Hy(li:ui-1,lj+1:uj,uk-1));

myxn(1,:,:) = -0.5* (DAT.Ez(li,lj:uj-1,lk:uk-1)+DAT.Ez(li,lj+1:uj,lk:uk-1));
mzxn(1,:,:) = 0.5 * (DAT.Ey(li,lj:uj-1,lk:uk-1)+DAT.Ey(li,lj:uj-1,lk+1:uk));
mxyn(:,1,:) = 0.5 * (DAT.Ez(li:ui-1,lj,lk:uk-1)+DAT.Ez(li+1:ui,lj,lk:uk-1));
mzyn(:,1,:) = -0.5* (DAT.Ex(li:ui-1,lj,lk:uk-1)+DAT.Ex(li:ui-1,lj,lk+1:uk));
mxzn(:,:,1) = -0.5* (DAT.Ey(li:ui-1,lj:uj-1,lk)+DAT.Ey(li+1:ui,lj:uj-1,lk));
myzn(:,:,1) = 0.5 * (DAT.Ex(li:ui-1,lj:uj-1,lk)+DAT.Ex(li:ui-1,lj+1:uj,lk));

jyxn(1,:,:) = 0.25 *(DAT.Hz(li,lj:uj-1,lk:uk-1)+DAT.Hz(li,lj:uj-1,lk+1:uk)+DAT.Hz(li-1,lj:uj-1,lk:uk-1)+DAT.Hz(li-1,lj:uj-1,lk+1:uk));
jzxn(1,:,:) = -0.25*(DAT.Hy(li,lj:uj-1,lk:uk-1)+DAT.Hy(li,lj+1:uj,lk:uk-1)+DAT.Hz(li-1,lj:uj-1,lk:uk-1)+DAT.Hz(li-1,lj+1:uj,lk:uk-1));
jzyn(:,1,:) = 0.25 *(DAT.Hx(li:ui-1,lj,lk:uk-1)+DAT.Hx(li+1:ui,lj,lk:uk-1)+DAT.Hx(li:ui-1,lj-1,lk:uk-1)+DAT.Hx(li+1:ui,lj-1,lk:uk-1));
jxyn(:,1,:) = -0.25*(DAT.Hz(li:ui-1,lj,lk:uk-1)+DAT.Hz(li:ui-1,lj,lk+1:uk)+DAT.Hz(li:ui-1,lj-1,lk:uk-1)+DAT.Hz(li:ui-1,lj-1,lk+1:uk));
jyzn(:,:,1) = -0.25*(DAT.Hx(li:ui-1,lj:uj-1,lk)+DAT.Hx(li+1:ui,lj:uj-1,lk)+DAT.Hx(li:ui-1,lj:uj-1,lk-1)+DAT.Hx(li+1:ui,lj:uj-1,lk-1));
jxzn(:,:,1) = 0.25 *(DAT.Hy(li:ui-1,lj:uj-1,lk)+DAT.Hy(li:ui-1,lj+1:uj,lk)+DAT.Hy(li:ui-1,lj:uj-1,lk-1)+DAT.Hy(li:ui-1,lj+1:uj,lk-1));

