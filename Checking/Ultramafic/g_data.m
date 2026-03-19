T_tc = [827.328 740.722] + 273.15;
P_tc = [25.4000e8   40e8];
g_tc = [-6399.52636 ,-6130.85286];
% site = {  xMgM1     xFeM1};
z_tc = [0.83209   0.16791
    0.79807   0.20193];
% ideal       gamma    activity        prop               µ0         RT ln a
grt_tc(:,:,1) =    [
   0.57611   1.0204    0.58784    0.83209  -6545.07373     -4.86126
   0.0047341   1.6403  0.0077654    0.16791  -5609.72113    -44.45049];
grt_tc(:,:,2) =    [
   0.50830   1.0395    0.52836    0.79807  -6310.50308     -5.37805
   0.0082341   1.8302   0.015070    0.20193  -5364.22682    -35.36303
    ];
p_tc     = squeeze(grt_tc(:,4,:))';
mu0_tc   = squeeze(grt_tc(:,5,:))';
a_tc     = squeeze(grt_tc(:,3,:))';
a_id_tc  = squeeze(grt_tc(:,1,:))';
gam_tc   = squeeze(grt_tc(:,2,:))';
RTlna_tc = squeeze(grt_tc(:,6,:))';
% a_id_tc.*gam_tc - a_tc

% 8.3144*T_tc(1)*log(a_id_tc.*gam_tc) - RTlna_tc*1e3
