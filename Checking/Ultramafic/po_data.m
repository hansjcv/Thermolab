T_tc = [389.197 532.875] + 273.15;
P_tc = [1e8   40e8];
g_tc = [-149.45391 ,-92.24613];
% site = { xfeM2      xVM2};
z_tc = [0.98659   0.01341
       0.91841   0.08159];
% ideal       gamma    activity        prop               µ0         RT ln a
po_tc(:,:,1) =    [
  0.84034  0.63027    0.52964    0.10732   -141.69348     -3.50003
  0.98659  0.99335    0.98003    0.89268   -149.85497     -0.11111];
po_tc(:,:,2) =    [
    0.98912  0.94421    0.93394    0.65273    -89.31632     -0.45798
     0.91841  0.81644    0.74982    0.34727    -94.96273     -1.92954];
p_tc     = squeeze(po_tc(:,4,:))';
mu0_tc   = squeeze(po_tc(:,5,:))';
a_tc     = squeeze(po_tc(:,3,:))';
a_id_tc  = squeeze(po_tc(:,1,:))';
gam_tc   = squeeze(po_tc(:,2,:))';
RTlna_tc = squeeze(po_tc(:,6,:))';
% a_id_tc.*gam_tc - a_tc

% 8.3144*T_tc(1)*log(a_id_tc.*gam_tc) - RTlna_tc*1e3
