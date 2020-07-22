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
ppm = 3.5;

% flag to save figures as pdfs
saveFig = 0;

if saveFig
    mkdir('figure_output');
end

%% simulate intra and extracellular signals

% tissue parameters
% intraC, extraC, tumorI, tumorE
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

% sequence parameters
offsetArray = [-5:0.1:5];
midPoint = (length(offsetArray)+1)/2;
ppm_MTRasym = offsetArray(midPoint:1:end);

% simulate the magnetization
PCarray = [20 1 20 1]; % phosphate concentration
M_P = calcM(pHarray, T1array, T2array, PCarray, ...
    AmineCarray, AmideCarray, CreatCarray, GlucoCarray, offsetArray, FA);
MTRasym_P = (M_P(midPoint:-1:1,:) - M_P(midPoint:1:end,:))*100;

PCarray = [10 10 10 10];
M_noP = calcM(pHarray, T1array, T2array, PCarray, ...
    AmineCarray, AmideCarray, CreatCarray, GlucoCarray, offsetArray, FA);
MTRasym_noP = (M_noP(midPoint:-1:1,:) - M_noP(midPoint:1:end,:))*100;

%% plot

ind = 1:2:length(offsetArray); % avoid too much data points on teh plot
% plot z pectrum
figure('rend','painters','pos',[10 10 300 300]);
plot(offsetArray(ind), M_P(ind,1),'o-','Color',colorArray{1}); hold on;
plot(offsetArray(ind), M_P(ind,2),'o-','Color',colorArray{3});
plot(offsetArray(ind), M_P(ind,3),'v-','Color',colorArray{1});
plot(offsetArray(ind), M_P(ind,4),'v-','Color',colorArray{3});
ylim([0 1.1]);
xticks(-5:1:5);

xlabel('Offset frequency (ppm)');
ylabel('Normalized magnetization');
% legend({'Normal - intracellular','Normal - extracellular',...
%     'Tumor - intracellular','Tumor - extracellular'},'box','off');
legend({'Normal(i)','Normal(e)',...
    'Tumor(i)','Tumor(e)'},'box','off','location','best');

if saveFig
    export_fig([exportFig_path '/z_spectrum.pdf'],'-m2','-transparent'); close;
end

%% plot MTRasym
ind = 1:2:length(ppm_MTRasym);
figure('rend','painters','pos',[10 10 300 300]);
plot(ppm_MTRasym(ind), MTRasym_P(ind,1),'o-','Color',colorArray{1}); hold on;
plot(ppm_MTRasym(ind), MTRasym_P(ind,2),'o-','Color',colorArray{3});
plot(ppm_MTRasym(ind), MTRasym_P(ind,3),'v-','Color',colorArray{1});
plot(ppm_MTRasym(ind), MTRasym_P(ind,4),'v-','Color',colorArray{3});
plot(ppm_MTRasym(ind), MTRasym_noP(ind,1),':','Color',colorArray{1}); hold on;
plot(ppm_MTRasym(ind), MTRasym_noP(ind,2),':','Color',colorArray{3});
plot(ppm_MTRasym(ind), MTRasym_noP(ind,3),':','Color',colorArray{1});
plot(ppm_MTRasym(ind), MTRasym_noP(ind,4),':','Color',colorArray{3});
ylim([0 3.5]); xlim([0 5]);
xticks(0:1:5);
line([3 3],[0 5],'Color',[0.7 0.7 0.7],'LineStyle','--');
line([3.5 3.5],[0 5],'Color',[0.7 0.7 0.7],'LineStyle','--');
line([1.9 1.9],[0 5],'Color',[0.7 0.7 0.7],'LineStyle','--');
line([1.28 1.28],[0 5],'Color',[0.7 0.7 0.7],'LineStyle','--');

xlabel('Offset frequency (ppm)');
ylabel('MTRasym (%)');
legend({'Normal(i)','Normal(e)',...
    'Tumor(i)','Tumor(e)','10 mM [PO_4]'},'box','off');

if saveFig
    export_fig([exportFig_path '/MTRasym_spectrum.pdf'],'-m2','-transparent'); close;
end

%% stacked bar plot
ratioMat = [3/4 1/2; 1/4 1/2]; % [0.76 0.52; 0.24 0.48]; from literature
% NormalI TumorI
% NormalE TumorE

MTRasym3_P = reshape(MTRasym_P(ppm_MTRasym == ppm,:),[2 2]).*ratioMat;
MTRasym3_noP = reshape(MTRasym_noP(ppm_MTRasym == ppm,:),[2 2]).*ratioMat;

figure('rend','painters','pos',[10 10 300 300]);
b = bar([MTRasym3_P]', 'stacked');
xticklabels({'Normal','Tumor','Normal','Tumor'});
xtickangle(0);
ylabel(['MTRasym at ' num2str(ppm,'%.1f') ' ppm(%)']); ylim([0 2]);
b(1).FaceColor = 'flat';
b(1).CData(:,:) = [repmat(colorArray{1},[2 1])];
b(1).EdgeColor = [1 1 1];
b(2).FaceColor = 'flat';
b(2).CData(:,:) = [repmat(colorArray{3},[2 1])];
b(2).EdgeColor = [1 1 1];
legend({'Intracellular','Extracellular'},'box','off');

% diff in contrast (I E)
diffMTR_P = sum(MTRasym3_P(:,2)) - sum(MTRasym3_P(:,1));
diffMTR_P_sep = MTRasym3_P(:,2) - MTRasym3_P(:,1);
fprintf(['diffMTR wP: %.4f Intracellular %.4f Extracellular %.4f '...
    'percentage of extracellular %.4f \n'],...
    diffMTR_P,diffMTR_P_sep(1),diffMTR_P_sep(2),diffMTR_P_sep(2)/diffMTR_P);

if saveFig
    export_fig([exportFig_path '/stack_barplot.pdf'],'-m2','-transparent'); close;
end

%% functions
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
    SystemParam.PoolConc   = [AmineCarray(tInd)  ...
        AmideCarray(tInd)  CreatCarray(tInd)    GlucoCarray(tInd)]; % mM
    SystemParam.PoolT1     = [1       1       1       1]; % s
    SystemParam.PoolT2     = [0.01    0.01    0.01    0.01]; % s
    
    SystemParam.WaterParam.T1 = T1array(tInd); % s
    SystemParam.WaterParam.T2 = T2array(tInd); % s
    
    Kamine = 58.31 + 3.45*10^10 * 10.^(pHarray(tInd)-13.6) + ...
        2.26*10^5 * APO4(tInd) + 7550 * BPO4(tInd);
    Kamide = 27.41 + 7.53*10^7 * 10.^(pHarray(tInd)-13.6) + ...
        8.52 * APO4(tInd) + 13.92 * BPO4(tInd);
    Kcreatine = 14.14 + 1.34*10^9 * 10.^(pHarray(tInd)-13.6) + ...
        630 * APO4(tInd) + 0 * BPO4(tInd);
    Kglucose = 3.93*10^9 * 10.^(pHarray(tInd)-13.6) + ...
        6.46*10^4 * APO4(tInd) + 2.35*10^4 * BPO4(tInd) + ...
        2.64*10^8 * 10.^(-pHarray(tInd));
    PoolKex    = [Kamine    Kamide    Kcreatine    Kglucose]; % Hz
    
    PulseParam.FA = 0;
    fprintf('M0 \n');
    M0 = CESTsim_function_5pool_wMTC...
        (SystemParam, PulseParam, 200, PoolKex);
    
    PulseParam.FA = FA;
    fprintf('offset ');
    for offInd = 1:length(offsetArray)
        fprintf('%.2f ',offsetArray(offInd));
        M(offInd,tInd) = ...
            CESTsim_function_5pool_wMTC...
            (SystemParam, PulseParam, offsetArray(offInd), PoolKex)/M0;
    end
    fprintf('\n');
    
end

end