clc; close; clear;

colorArray = {[0, 0.4470, 0.7410],...
    [0.8500, 0.3250, 0.0980],...
    [0.9290, 0.6940, 0.1250],...
    [0.4940, 0.1840, 0.5560]};

exportFig_path = 'figure_output';
reset(0);
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultLineLineWidth',2);

%% !!! parameter to change for different CEST contrast !!!

% saturation pulse flip angle
FA = 1030;
% FA 1030 = peak B1 1.402uT
% FA 3090 = peak B1 4.206uT

% offset of interest
ppm = 1.2;

% flag to save figures as pdfs
saveFig = 0;

if saveFig
    mkdir('figure_output');
end

%% simulate dynamic function

timestep = 2; % min
T = 60; % min
t = 0:timestep:T; % time array

% Cp from literature 
% Dynamic imaging and tracer kinetic modeling for emission tomography using rotating detectors
A1 = 851.1; A2 = 21.88; A3 = 20.81; % uCi/ml
lambda1 = -4.134; lambda2 = -0.1191; lambda3 = -0.0104; % /min
Cp = (A1*t-A2-A3).*exp(lambda1*t) + A2*exp(lambda2*t) + A3*exp(lambda3*t);

% kinetics parameters from literature
% Noninvasive determination of local cerebral metabolic rate of glucose in man
% k = [0.0540 0.1090 0.0450 0.0058 0.0157]; % WM parameters
% k = [0.1020 0.1300 0.0620 0.0068 0.0329]; % GM parameters

% normal brain tissue
kBrain = [0.0540 0.1090 0.05 0.005 0.0450];
[~,y] = ode45(@(a,y) odefcn(a,y,Cp,t,kBrain),t,[0; 0]);
CfWM = y(:,1); 
CmWM = y(:,2);

figure('rend','painters','pos',[10 10 300 300]);
plot(t,Cp,'Color',colorArray{4}); hold on; % blood concentration
plot(t,CfWM,'Color',colorArray{3}); % interstitial concentration
plot(t,CmWM,'Color',colorArray{1}); % intracellular concentration
plot(t,CfWM+CmWM,'Color',colorArray{2}); % total tissue concentration

% tumor tissue
kTumor = [0.1020 0.1300 0.1 0.01 0.0620];
[~,y] = ode45(@(a,y) odefcn(a,y,Cp,t,kTumor),t,[0; 0]);
CfGM = y(:,1); CmGM = y(:,2);

plot(t,CfGM,':','Color',colorArray{3});
plot(t,CmGM,':','Color',colorArray{1});
plot(t,CfGM+CmGM,':','Color',colorArray{2});
ylabel('Concentration (A.U.)');
xlabel('Time (min)');
legend({'Arterial input function','Extracellular (Normal)',...
    'Intracellular (Normal)','Total tissue (Normal)','Tumor'},'box','off');

if saveFig
    export_fig([exportFig_path '/DGE_kinetics.pdf'],'-m2','-transparent'); close;
end

%% set up parameters

% tissue parameters
% normalI, normalE, tumorI, tumorE
pHarray = [7.0 7.4 7.2 6.6];
% normal WM https://www.ncbi.nlm.nih.gov/pubmed/10232510
% glioma T2 https://onlinelibrary.wiley.com/doi/pdf/10.1002/jmri.20335
% glioma T1 http://www.ajnr.org/content/ajnr/early/2016/12/29/ajnr.A5035.full.pdf
T1array = [0.832 0.832 1.558 1.558]; 
T2array = [0.110 0.110 0.160 0.160];
AmineCarray = [20 10 24 20];
AmideCarray = [70 15 84 18];
CreatCarray = [20 0.1 16 0.1];
GlucoCarray = [15 15 30 30];
GlucoCarray_external = [CmWM CfWM CmGM CfGM];

%% Simulate DGE signal

% S0
PCarray = [20 1 20 1]; % phosphate concentrations
M0_P = calcM(pHarray, T1array, T2array, PCarray, ...
        AmineCarray, AmideCarray, CreatCarray, ...
        GlucoCarray, ...
        0, 0);
    
% S(t)
M_P = zeros(size(GlucoCarray_external));
for i = 1:length(t)
    fprintf('%d \n',i);
    M_P(i,:) = calcM(pHarray, T1array, T2array, PCarray, ...
        AmineCarray, AmideCarray, CreatCarray, ...
        GlucoCarray+GlucoCarray_external(i,:), ...
        ppm, FA);
end

S_P = M_P./M0_P*100;

%% plot deltaS/s0
deltaM_P = (repmat(M_P(1,:),length(t),1) - M_P)./M0_P*100;

figure('rend','painters','pos',[10 10 300 300]);
plot(t, deltaM_P(:,1),'o-','Color',colorArray{1}); hold on;
plot(t, deltaM_P(:,2),'o-','Color',colorArray{3});
plot(t, deltaM_P(:,1)+deltaM_P(:,2),'o-','Color',colorArray{2});
plot(t, deltaM_P(:,3),'v-','Color',colorArray{1});
plot(t, deltaM_P(:,4),'v-','Color',colorArray{3});
plot(t, deltaM_P(:,3)+deltaM_P(:,4),'v-','Color',colorArray{2});

xlabel('Time (min)');
ylabel('\DeltaS/S_0 (%)'); ylim([0 0.85]);
legend({'Normal(i)','Normal(e)','Normal(total)',...
    'Tumor(i)','Tumor(e)','Tumor(total)'},'box','off');

if saveFig
    export_fig([exportFig_path '/DGE_signals.pdf'],'-m2','-transparent'); close;
end

%% functions
function [dydt] = odefcn(a,y,Cp,t,k)

Cp_value = interp1(t,Cp,a);
dydt = zeros(2,1);
dydt(1) = k(1)*Cp_value-(k(2)+k(3))*y(1)+k(4)*y(2);
dydt(2) = k(3)*y(1)-(k(4)+k(5))*y(2);

end

function [M] = calcM(pHarray, T1array, T2array, PCarray, ...
    AmineCarray, AmideCarray, CreatCarray, GlucoCarray, offsetArray, FA)

pKa = 6.82;
Kratio = 10.^(pHarray-pKa);
BPO4 = PCarray./(1+Kratio)/1000;
APO4 = PCarray/1000 - BPO4;

SystemParam.MTCParam.Cmacro = 0;  % percentage
SystemParam.MTCParam.Kmacro = 20;  % Hz
SystemParam.MTCParam.T1 = 1;           % s
SystemParam.MTCParam.T2 = 10*10^-6;    % s

M = zeros([length(offsetArray) length(pHarray)]);
for tInd = 1:length(pHarray)
    
    SystemParam.PoolOffset = [3.0     3.5     1.9     1.3]; % ppm
    SystemParam.PoolConc   = [AmineCarray(tInd)  AmideCarray(tInd)...
                              CreatCarray(tInd)  GlucoCarray(tInd)]; % mM
    SystemParam.PoolT1     = [1       1       1       1]; % s
    SystemParam.PoolT2     = [0.01    0.01    0.01    0.01]; % s
    
    SystemParam.WaterParam.T1 = T1array(tInd); % s
    SystemParam.WaterParam.T2 = T2array(tInd); % s
    
    Kamine = 58.2 + 3.45*10^10 * 10.^(pHarray(tInd)-13.6) + ...
        2.26*10^5 * APO4(tInd) + 7550 * BPO4(tInd);
    Kamide = 27.4 + 7.53*10^7 * 10.^(pHarray(tInd)-13.6) + ...
        8.51 * APO4(tInd) + 13.92 * BPO4(tInd);
    Kcreatine = 14.14 + 1.34*10^9 * 10.^(pHarray(tInd)-13.6) + ...
        630 * APO4(tInd) + 0 * BPO4(tInd);
    Kglucose = 3.93*10^9 * 10.^(pHarray(tInd)-13.6) + ...
       6.46*10^4 * APO4(tInd) + 2.35*10^4 * BPO4(tInd) + ...
       1.64*10^8 * 10.^(-pHarray(tInd));
    PoolKex    = [Kamine    Kamide    Kcreatine    Kglucose]; % Hz
    
%     PulseParam.FA = 0;
%     fprintf('M0 \n');
%     M0 = CESTsim_function_5pool_wMTC...
%         (SystemParam, PulseParam, 200, PoolKex);
    
    PulseParam.FA = FA;
%     fprintf('offset ');
    for offInd = 1:length(offsetArray)
%         fprintf('%.2f ',offsetArray(offInd));
        M(offInd,tInd) = ...
            CESTsim_function_5pool_wMTC...
            (SystemParam, PulseParam, offsetArray(offInd), PoolKex);
    end
%     fprintf('\n');
    
end

end