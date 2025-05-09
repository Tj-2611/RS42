function [V0, V_bd, options, optimInfo, paramDefCell, ...
    objFctSettings, objFctParams, DEParams, emailParams] = initial_V(model, opt)
% Returns the options, initializations and ranges for variance component optimization. 

if strcmp(opt,'ip')
    options = optimoptions('fmincon','Algorithm','interior-point',...
'MaxFunctionEvaluations',10^4,'MaxIterations',10^4,'FunctionTolerance',10^(-12),...
'StepTolerance',10^(-4),'FiniteDifferenceStepSize',10^(-4),'Display','iter-detailed');
elseif strcmp(opt,'lm')
    options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt',...
    'MaxFunctionEvaluations',10^4,'MaxIterations',10^4,'FunctionTolerance',10^(-12),...
    'StepTolerance',10^(-4),'FiniteDifferenceStepSize',10^(-4),'Display','iter-detailed');
end

if strcmp(model,'H2Fwithmu')
    V0 = [3e-2; 3e-2];
    V_bd = [[1e-4, 1.5]; [1e-4, 1.5]];
    quantization = [1, 1]*1e-3;
elseif strcmp(model,'42withmu')     
    V0 = [3e-3; 3e-3];
    V_bd = [[0, 1]; [0, 1]];    
    quantization = [1, 1]*1e-3;    
elseif strcmp(model,'42rswithmu')      
    V0 = [3e-2; 3e-2];
    V_bd = [[1e-4, 1.5]; [1e-4, 1.5]];
    quantization = [1, 1]*1e-3;
end


optimInfo.title = [model,'_', opt];
paramDefCell = {
	'', V_bd, quantization, V0
};
objFctSettings = {};
objFctParams = [];
% Get default DE parameters
DEParams = getdefaultparams;
% Set number of population members (often 10*D is suggested) 
DEParams.NP = 10*length(V0);
DEParams.feedSlaveProc = 0;

% Set times
DEParams.maxiter  = 5e2/2*length(V0);
DEParams.maxtime  = 1200; % in seconds
DEParams.maxclock = [];
% Set display options
DEParams.infoIterations = 10;
DEParams.infoPeriod     = 100; % in seconds
DEParams.saveHistory = false;
% Do not send E-mails
emailParams = [];
end

