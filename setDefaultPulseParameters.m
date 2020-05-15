function [PulseParam] = setDefaultPulseParameters(PulseParam)

% Pulse parameters --------------------------------------------------------
% - FieldStrength
% - BaseFrequency
% - B1peak, FA
% - GaussPulse
% - PulseLength
% - TimeStep
% - ReadoutLength
% - EchoTime
% - Deadtime
% - NumberOfPulses
% - Nslice
% -------------------------------------------------------------------------

GMR = 42.576;           % MHz/T

if ~isfield(PulseParam,'FieldStrength')
    PulseParam.FieldStrength = 3.0;    % T
end

PulseParam.BaseFrequency = GMR * PulseParam.FieldStrength;  % MHz

if ~isfield(PulseParam,'FA')
    PulseParam.FA = 3090;    % degree
end

if ~isfield(PulseParam,'PulseLength')
    PulseParam.PulseLength = 99840/1000000;    % s RF pulse length
end

GaussPulseNorm = [0.037775253	0.043008072	0.048837602	0.055312138	0.062481052	0.070394379	0.079102353	0.088654867	0.099100885	0.11048778...
    0.122860627	0.13626144	0.150728364	0.166294837	0.182988716	0.200831402	0.219836953	0.240011216	0.261350988	0.283843216	0.307464274	0.3321793...
    0.357941651	0.38469246	0.412360329	0.440861168	0.470098193	0.499962099	0.53033141	0.561073017	0.592042906	0.623087081	0.654042662	0.684739173...
    0.714999991	0.744643949	0.773487072	0.801344424	0.828032047	0.853368951	0.87717914	0.899293634	0.919552457	0.937806557	0.953919631	0.967769819...
    0.979251234	0.988275312	0.994771952	0.998690418	1	0.998690418	0.994771952	0.988275312	0.979251234	0.967769819	0.953919631	0.937806557...
    0.919552457	0.899293634	0.87717914	0.853368951	0.828032047	0.801344424	0.773487072	0.744643949	0.714999991	0.684739173	0.654042662	0.623087081...
    0.592042906	0.561073017	0.53033141	0.499962099	0.470098193	0.440861168	0.412360329	0.38469246	0.357941651	0.3321793	0.307464274	0.283843216...
    0.261350988	0.240011216	0.219836953	0.200831402	0.182988716	0.166294837	0.150728364	0.13626144	0.122860627	0.11048778	0.099100885	0.088654867...
    0.079102353	0.070394379	0.062481052	0.055312138	0.048837602	0.043008072	0.037775253];
PulseParam.PulseBins = length(GaussPulseNorm);
PulseParam.TimeStep = PulseParam.PulseLength / PulseParam.PulseBins;    

if ~isfield(PulseParam,'B1peak')
    PulseParam.B1peak = PulseParam.FA./sum(GaussPulseNorm*GMR*PulseParam.TimeStep*360);
end
PulseParam.GaussPulse = GaussPulseNorm .* (PulseParam.B1peak*2*pi*GMR);

if ~isfield(PulseParam,'ReadoutLength')
    PulseParam.ReadoutLength = 0.118;    % in seconds
end

if ~isfield(PulseParam,'EchoTime')
    PulseParam.EchoTime = 0.015;    % in seconds
end

if ~isfield(PulseParam,'Deadtime')
    PulseParam.Deadtime = 0.005;    % in seconds
end

if ~isfield(PulseParam,'NumberOfPulses')
    PulseParam.NumberOfPulses = [3 3]; % primary & secondary
end

if ~isfield(PulseParam,'Nslice')
    PulseParam.Nslice = 10; % number of slices (excitations) to loop through
end

end