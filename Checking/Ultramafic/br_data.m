T_tc = [389.197 532.875] + 273.15;
P_tc = [1e8   40e8];
g_tc = [-938.08270,-850.29871];
% site = { x(Mg)     x(Fe)};
z_tc = [0.86697   0.13303
        0.82140   0.17860];
% ideal       gamma    activity        prop               µ0         RT ln a
br_tc(:,:,1) =    [
    0.86697        0    0.86697    0.86697   -979.81880     -0.78610
    0.13303        0    0.13303    0.13303   -649.84165    -11.10883 ];
br_tc(:,:,2) =    [
    0.82140        0    0.82140    0.82140   -906.11062     -1.31854
    0.17860        0    0.17860    0.17860   -576.01293    -11.54409];
p_tc     = squeeze(br_tc(:,4,:))';
mu0_tc   = squeeze(br_tc(:,5,:))';
a_tc     = squeeze(br_tc(:,3,:))';
a_id_tc  = squeeze(br_tc(:,1,:))';
gam_tc   = squeeze(br_tc(:,2,:))';
RTlna_tc = squeeze(br_tc(:,6,:))';
% a_id_tc.*gam_tc - a_tc

% 8.3144*T_tc(1)*log(a_id_tc.*gam_tc) - RTlna_tc*1e3
