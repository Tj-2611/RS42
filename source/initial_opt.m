function [paramdl0, paramdl_bd, options, optimInfo, paramDefCell, ...
    objFctSettings, objFctParams, DEParams, emailParams] = initial_opt(model, opt)
% Returns the options, initializations and ranges for structural parameter optimization. 

if strcmp(opt,'ip')
    options = optimoptions('fmincon','Algorithm','interior-point',...
'MaxFunctionEvaluations',10^4,'MaxIterations',30,'FunctionTolerance',10^(-12),...
'StepTolerance',10^(-4),'FiniteDifferenceStepSize',10^(-4),'Display','iter-detailed');
elseif strcmp(opt,'lm')
    options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt',...
    'MaxFunctionEvaluations',10^4,'MaxIterations',30,'FunctionTolerance',10^(-12),...
    'StepTolerance',10^(-4),'FiniteDifferenceStepSize',10^(-4),'Display','iter-detailed');
end

if strcmp(model, 'H2Fwithmu')
     paramdl0 = [0.5, 0.2, 0.005, 2, -0.75, ...
                  0.2, 0.15, 0.2, -0.75]';
    paramdl_bd = [[0, 1]; [0, 50]; [0, 10]; [0, 50]; [-1, 1]; ...
                          [0, 50]; [0, 10]; [0, 50]; [-1, 1]]; 
    quantization = zeros(size(paramdl_bd, 1), 1);
elseif strcmp(model,'42withmu') 
    paramdl0 = [0.05, 1, 0.5, 0.1, -0.5, ...
                      50, 0.2, 70, -0.5]'; 
    paramdl_bd = [[0, 1]; [1e-4, 50]; [1e-4, 10]; [1e-4, 50]; [-1, 1]; ...
                          [1e-2, 100]; [1e-4, 10]; [10, 100]; [-1, 1]];
    quantization = zeros(size(paramdl_bd, 1), 1); 
elseif strcmp(model,'42rswithmu') 
    paramdl0 = [0.05, 0.1, 1, 0.5, 0.1, -0.5, ...
                      50, 0.2, 70, -0.5]'; 
    paramdl_bd = [[0, 1]; [1e-4, 0.4999];[1e-4, 50]; [1e-4, 10]; [1e-4, 50]; [-1, 1]; ...
                                         [1e-2, 100]; [1e-4, 10]; [10, 100]; [-1, 1]];
    quantization = zeros(size(paramdl_bd, 1), 1);
end


optimInfo.title = [model,'_', opt];
paramDefCell = {
	'', paramdl_bd, quantization, paramdl0
};
objFctSettings = {};
objFctParams = [];
% Get default DE parameters
DEParams = getdefaultparams;
% Set number of population members (often 10*D is suggested) 
DEParams.NP = 10*length(paramdl0);
DEParams.feedSlaveProc = 0;

% Set times
DEParams.maxiter  = 5e2/2*length(paramdl0);
DEParams.maxtime  = 3600/6*length(paramdl0); % in seconds
DEParams.maxclock = [];
% Set display options
DEParams.infoIterations = 10;
DEParams.infoPeriod     = 100; % in seconds
DEParams.saveHistory = false;
% Do not send E-mails
emailParams = [];

end
