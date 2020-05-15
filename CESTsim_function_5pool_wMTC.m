function [M_sim, PulseParam] = ...
    CESTsim_function_5pool_wMTC(SystemParam, PulseParam, OffsetFreq, PoolKex)

% kex = k0 + kb * 10^(pH-14) + kcB[H2PO4 -] + kcA[HPO4 2-]

% Input -------------------------------------------------------------------
% SystemParam
%   PoolOffset A, B, C, D
%   PoolConc
%   PoolT1
%   PoolT2
%   PoolKex
%   WaterParam
%       T1, T2
%   MTCParam
%       Cmacro, Kmacro, T1, T2
% PulseParam

[PulseParam] = setDefaultPulseParameters(PulseParam);

%% System parameters

if nargin > 3
    SystemParam.PoolKex = PoolKex;
end

[~, M0List, KxwList, KwxList, T1List, T2List, ~, ~, dList] ...
    = unpackSystemParameters(SystemParam, PulseParam, OffsetFreq);

%% Pre-calculations

% Initial magnetization
x0 = [0;0;M0List(1);0;0;M0List(2);0;0;M0List(3);0;0;M0List(4);0;0;M0List(5);M0List(6)]; 

% MTC component
if M0List(6) > 1e-15
    fun  = @(x) ...
        sin(x).*sqrt(2/pi)*T2List(6)./(abs(3*cos(x).^2-1)).*exp(-2*(dList(6)*T2List(6)./(abs(3*cos(x).^2-1))).^2);
    g = integral(fun,0,pi/2);
else
    g = 0;
end

[b, B, Binvb] = ...
    calculateBlochMatrixB(M0List, KxwList, KwxList, T1List, T2List, dList);

expB_TimeStep = expm(B*PulseParam.TimeStep);
expB_ReadoutLength = expm(B*PulseParam.ReadoutLength);
expB_EchoTime = expm(B*PulseParam.EchoTime);
expB_DeadTime = expm(B*PulseParam.Deadtime);

expAList = zeros([PulseParam.PulseBins size(expB_TimeStep) ]);
AinvbList = zeros([PulseParam.PulseBins size(Binvb)]);
for z = 1:PulseParam.PulseBins  % z is Gaussian Pulse shape loop
    % Input variable B1 amplitude (use longGaussPulse if PulseBin ~= linspace desired)
    w_1 = PulseParam.GaussPulse(z);
    [~, A, Ainvb] = ...
        calculateBlochMatrixA(M0List, KxwList, KwxList, T1List, T2List, dList, w_1, g);
    expAList(z,:,:) = expm(A*PulseParam.TimeStep);
    AinvbList(z,:,:) = Ainvb;
end

%% Loop through time steps

for m = 1:PulseParam.Nslice  % m is multiple RF influence from short TR loop (#slices)
    
    if m == 1
        NumberOfPulses = PulseParam.NumberOfPulses(1);
    else
        NumberOfPulses = PulseParam.NumberOfPulses(2);
    end
    
    % RF saturation
    for N_RF = 1:NumberOfPulses  %% NumberOfPulses is number of RF pulses in the pulse train
        
        for z = 1:PulseParam.PulseBins  % z is Gaussian Pulse shape loop
            
            % FROM SHERRY CODE doi:  10.1002/nbm.3237
            if PulseParam.GaussPulse(z) == 0
                y0 = x0;
                x0 = expB_TimeStep*(y0+Binvb)-Binvb;
            else
                y0 = x0;
                x0 = squeeze(expAList(z,:,:))*(y0+AinvbList(z,:)')-AinvbList(z,:)';
            end
            
        end
        
        y0 = x0;
        x0 = expB_DeadTime*(y0+Binvb)-Binvb;
        
    end
    
    % Spoiling
    x0([1 2 4 5 7 8 10 11 13 14]) = 0;
    
    % Water only excitation (assume magnetization flip into xy plane: Mx)
    if m == PulseParam.Nslice
        x0(1) = x0(3); x0(3) = 0;
        
        % This accounts for decay during TE
        y0 = x0;
        x0 = expB_EchoTime*(y0+Binvb)-Binvb;
        M_sim = sqrt(x0(1)^2+x0(2)^2);
%         M_sim = x0(1);
        
    end
    
    % This accounts for decay during readout
    y0 = x0;
    x0 = expB_ReadoutLength*(y0+Binvb)-Binvb;
    %         x0(7) = M_a0 - ( M_a0 - x0(7) ) * exp( - ReadoutLength / T_1a );
    %         x0(8) = M_b0 - ( M_b0 - x0(8) ) * exp( - ReadoutLength / T_1b );
    %         x0(9) = M_c0 - ( M_c0 - x0(9) ) * exp( - ReadoutLength / T_1c );
    
end

