% Calculate total radiated power
radiated_power = 0;

powr = dx*dy*sum(myzp.*conj(jxzp)-mxzp.*conj(jyzp),"all");
powr = powr-dx*dy*sum(myzn.*conj(jxzn)-mxzn.*conj(jyzn),"all");
powr = powr+dx*dz*sum(mxyp.*conj(jzyp)-mzyp.*conj(jxyp),"all");
powr = powr-dx*dz*sum(mxyn.*conj(jzyn)-mzyn.*conj(jxyn),"all");
powr = powr+dy*dz*sum(mzxp.*conj(jyxp)-myxp.*conj(jzxp),"all");
powr = powr-dy*dz*sum(mzxn.*conj(jyxn)-myxn.*conj(jzxn),"all");
radiated_power = 0.5*real(powr);