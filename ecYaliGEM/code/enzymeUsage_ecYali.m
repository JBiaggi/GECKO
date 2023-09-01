adapterLocation = fullfile(findGECKOroot,'ecYaliGEM','ecYaliGEMAdapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();
ecModel = loadEcModel('ecYaliGEM_stage4.yml');

sol = solveLP(ecModel,1);

usageData = enzymeUsage(ecModel,sol.x);

usageReport = reportEnzymeUsage(ecModel, usageData);

%ecModel = loadEcModel('ecYaliGEM_FSEOF.yml');
ecModel = loadEcModel('ecYaliGEM_FSEOF_pooled.yml');

%% from ecFactory
CS_index = find(strcmpi(ecModel.rxns,'y001808'));
growthPos = find(strcmpi(ecModel.rxns,'xBIOMASS'));
CS_MW = 0.09209;

%Enable cellular growth
ecModel = setParam(ecModel,'lb',growthPos, 0);
ecModel = setParam(ecModel,'ub',growthPos, 1000);

%Set a fixed unit glucose uptake rate
ecModel = setParam(ecModel,'ub',CS_index,1);

%Get biomass yield for a unit glycerol uptake rate
ecModel = setParam(ecModel,'obj',growthPos,1);
solution      = solveLP(ecModel,1);
WT_yield      = -solution.x(growthPos)/(solution.x(CS_index)*CS_MW);
disp(['* The maximum biomass yield is ' num2str(WT_yield) '[g biomass/g carbon source]']);

%Obtain a suboptimal yield value to run ecFactory
expYield = 0.49*WT_yield;
disp('* The ecFactory method will scan flux distributions spanning from')
disp(['a suboptimal biomass yield of: ' num2str(0.5*expYield) ' to: ' num2str(2*expYield) ' [g biomass/g carbon source]']);

%%
FC = ecFSEOF(ecModel,'EXC_OUT_m1640','y001808',[0.5*expYield 2*expYield],[],[]);

%FC = ecFSEOF(ecModel,'EXC_OUT_m1727','y001808',[0.1 0.9]);



