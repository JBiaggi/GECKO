adapterLocation = fullfile(findGECKOroot,'ecYali4GEM','ecYali4Adapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();
ecModel = loadEcModel('ecYali4GEM.yml');

sol = solveLP(ecModel);
ecModel = setParam(ecModel,'eq','1808',sol.x(find(strcmpi(ecModel.rxns,'1808'))));
ecModel = setParam(ecModel,'eq','biomass_C',-sol.f*0.90);
ecModel = setParam(ecModel, 'obj', 'prot_pool_exchange', 1);
sol   = solveLP(ecModel);
ecModel = setParam(ecModel, 'lb', 'prot_pool_exchange', sol.x(strcmpi(ecModel.rxns, 'prot_pool_exchange')) * 1.01);
ecModel = setParam(ecModel,'obj','EXC_OUT_m1640',1);
sol = solveLP(ecModel,1);

%printFluxes(ecModel,sol.x,false,[],'fluxes_ecYali_SBY145_Nlim_MaxLipids.json','"%rxnID": %flux,')

usageData = enzymeUsage(ecModel,sol.x);

usageReport = reportEnzymeUsage(ecModel, usageData);

%% from ecFactory
ecModel = loadEcModel('ecYali4GEM.yml');
%ecModel = loadEcModel('ecYaliGEM_SBY145_Nlim.yml');

CS_index = find(strcmpi(ecModel.rxns,'1808'));
growthPos = find(strcmpi(ecModel.rxns,'biomass_C'));
CS_MW = 0.09209;

%Enable cellular growth
ecModel = setParam(ecModel,'lb',growthPos, 0);
ecModel = setParam(ecModel,'ub',growthPos, 1000);

%Set a fixed unit glycerol uptake rate
ecModel = setParam(ecModel,'eq',CS_index,-1);

%Get biomass yield for a unit glycerol uptake rate
ecModel = setParam(ecModel,'obj',growthPos,1);
solution      = solveLP(ecModel,1);
WT_yield      = -solution.x(growthPos)/(solution.x(CS_index)*CS_MW);
disp(['* The maximum biomass yield is ' num2str(WT_yield) '[g biomass/g carbon source]']);

%Obtain a suboptimal yield value to run ecFactory
expYield = 0.463;
disp('* The ecFactory method will scan flux distributions spanning from')
disp(['a suboptimal biomass yield of: ' num2str(0.5*expYield) ' to: ' num2str(0.9*WT_yield) ' [g biomass/g carbon source]']); %ecFactory
%disp(['a suboptimal biomass yield of: ' num2str(expYield) ' to: ' num2str(0.9*WT_yield) ' [g biomass/g carbon source]']); %JSB

%%
%FC = ecFSEOF(ecModel,'EXC_OUT_m1640','1808',[0.5*expYield/WT_yield 0.9],[],[]); %ecFactory
FC = ecFSEOF(ecModel,'EXC_OUT_m1640','1808');
%FC = ecFSEOF(ecModel,'EXC_OUT_m1640','y001808',[0.736825014 0.74139142],[],[]); %JSB



