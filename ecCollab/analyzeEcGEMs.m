adapterLocation = fullfile(findGECKOroot,'ecCollab','ecCollabAdapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();


% Load model: already constrained with 10% variance around chemostat fluxes
ecModel = loadEcModel('ecJFYL07_prot.yml');

% Set bounds for biomass too
sol = solveLP(ecModel);
ecModel = setParam(ecModel,'var', 'xBIOMASS', -sol.f, 10);

% Get good reactions
[~, goodRxns] = randomSampling(ecModel,1,true,true,true);

% Run random sampling
solutions = randomSampling(ecModel,10000,true,true,true,goodRxns);

% Data treatment
fluxMean = full(mean(solutions,2));
fluxSD = full(std(solutions,0,2));
logVec = abs(fluxMean) > abs(fluxSD);
fluxMeanSD = logVec.*fluxMean;

% Export data
fluxesForEscher(ecModel.rxns,fluxMeanSD,'ecJFYL07_prot_FBA_SD.json');