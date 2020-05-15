function [wList, M0List, KxwList, KwxList, T1List, T2List, k1List, k2List, dList] ...
    = unpackSystemParameters(SystemParam, PulseParam, OffsetFreq)

% exchange pool offsets
w_w = 0.0; % PulseParam.BaseFrequency * SystemParam.PoolOffset(1);  % bulk water frequency
w_a = PulseParam.BaseFrequency * SystemParam.PoolOffset(1);  % Hz 
w_b = PulseParam.BaseFrequency * SystemParam.PoolOffset(2);  % Hz 
w_c = PulseParam.BaseFrequency * SystemParam.PoolOffset(3);  % Hz 
w_d = PulseParam.BaseFrequency * SystemParam.PoolOffset(4);  % Hz 
w_m = 0.0; % PulseParam.BaseFrequency * SystemParam.PoolOffset(6);  % Hz 
wList = [w_w w_a w_b w_c w_d w_m];
% w: water pool, m: semi-solid pool

% exchange pool concentrations
M_a0 = SystemParam.PoolConc(1)/(55.5*2)/1000; % mM 
M_b0 = SystemParam.PoolConc(2)/(55.5*2)/1000; % mM 
M_c0 = SystemParam.PoolConc(3)/(55.5*2)/1000; % mM 
M_d0 = SystemParam.PoolConc(4)/(55.5*2)/1000; % mM 
M_m0 = SystemParam.MTCParam.Cmacro; % percentage
M_w0 = 1; % - M_a0 - M_b0 - M_c0 - M_d0; % water 55.5M
M0List = [M_w0 M_a0 M_b0 M_c0 M_d0 M_m0];
M0List(M0List < 1e-16) = 0;

% effective exchange rate
k_aw = SystemParam.PoolKex(1); % Hz 
k_bw = SystemParam.PoolKex(2); % Hz 
k_cw = SystemParam.PoolKex(3); % Hz 
k_dw = SystemParam.PoolKex(4); % Hz 
k_mw = SystemParam.MTCParam.Kmacro; % Hz 

k_wa = (M_a0/M_w0)*k_aw;
k_wb = (M_b0/M_w0)*k_bw;
k_wc = (M_c0/M_w0)*k_cw;
k_wd = (M_d0/M_w0)*k_dw;
k_wm = (M_m0/M_w0)*k_mw;

k_w = k_wa + k_wb + k_wc + k_wd + k_wm;
KxwList = [k_w k_aw k_bw k_cw k_dw k_mw];
KwxList = [k_w k_wa k_wb k_wc k_wd k_wm];

% relaxation rates
T_1w = SystemParam.WaterParam.T1; % 3.374 from Chen et al % 1.22 for tissue % 1.36 for tumor (Chistopher Larsson 2015) %% 1.939 7T
T_2w = SystemParam.WaterParam.T2; % 2.5 from Araujo et al % 0.107 for tissue (Ellingson paper) % 0.170 for NET (Ellingson paper) % 0.055 7T

T_1a = SystemParam.PoolT1(1);
T_1b = SystemParam.PoolT1(2);
T_1c = SystemParam.PoolT1(3);
T_1d = SystemParam.PoolT1(4);
T_1m = SystemParam.MTCParam.T1; % 1s
T1List = [T_1w T_1a T_1b T_1c T_1d T_1m];

T_2a = SystemParam.PoolT2(1);
T_2b = SystemParam.PoolT2(2);
T_2c = SystemParam.PoolT2(3);
T_2d = SystemParam.PoolT2(4);
T_2m = SystemParam.MTCParam.T2; % 10*10^-6s
T2List = [T_2w T_2a T_2b T_2c T_2d T_2m];

k_1w = (1 / T_1w) + k_w;
k_2w = (1 / T_2w) + k_w;
k_1a = (1 / T_1a) + k_aw;
k_2a = (1 / T_2a) + k_aw;
k_1b = (1 / T_1b) + k_bw;
k_2b = (1 / T_2b) + k_bw;
k_1c = (1 / T_1c) + k_cw;
k_2c = (1 / T_2c) + k_cw;
k_1d = (1 / T_1d) + k_dw;
k_2d = (1 / T_2d) + k_dw;
k_1m = (1 / T_1m) + k_mw;
k_2m = (1 / T_2m) + k_mw;

k1List = [k_1w k_1a k_1b k_1c k_1d k_1m];
k2List = [k_2w k_2a k_2b k_2c k_2d k_2m];

% saturation offset
w = PulseParam.BaseFrequency * OffsetFreq;  % applied RF saturation frequency
dw = (w - w_w) * 2*pi;
da = (w - w_a) * 2*pi;
db = (w - w_b) * 2*pi;
dc = (w - w_c) * 2*pi;
dd = (w - w_d) * 2*pi;
dm = (w - w_m) * 2*pi;
dList = [dw da db dc dd dm];

% remove pools not used
M0List(M0List == 0) = 0;
KxwList(M0List == 0) = 0;
KwxList(M0List == 0) = 0;

end