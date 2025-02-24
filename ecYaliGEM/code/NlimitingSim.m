adapterLocation = fullfile(findGECKOroot,'ecYaliGEM','ecYaliGEMAdapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();

CNratios = [10 20 30 40 50 60];

ecModel = loadEcModel('ecYaliGEM.yml');

CS_index = find(strcmpi(ecModel.rxns,'y001808')); % glycerol
growthPos = find(strcmpi(ecModel.rxns,'xBIOMASS'));
NsourcePos = find(strcmpi(ecModel.rxns,'y001654')); % NH4

% Adjust glycerol uptake to N-limiting phase compatible value
ecModel = setParam(ecModel,'eq','y001808',-0.47);

% Set objective for lipids production with constrained biomass
sol = solveLP(ecModel);
ecModel = setParam(ecModel,'eq','xBIOMASS',-sol.f*0.99);
ecModel = setParam(ecModel, 'obj', 'prot_pool_exchange', 1);
sol   = solveLP(ecModel);

printFluxes(ecModel,sol.x)

ecModel = setParam(ecModel, 'lb', 'prot_pool_exchange', sol.x(strcmpi(ecModel.rxns, 'prot_pool_exchange')) * 1.01);
ecModel = setParam(ecModel,'obj','EXC_OUT_m1640',1);
sol = solveLP(ecModel);

%printFluxes(ecModel,sol.x)

CNratios = [sol.x(CS_index)*3/sol.x(NsourcePos) CNratios];

growthRates = zeros(1,length(CNratios));
growthRates(1) = sol.x(growthPos);

lipidProduction = zeros(1,length(CNratios));
lipidProduction(1) = -sol.f;

for i = 2:length(CNratios)
    ecModel = setParam(ecModel,'lb',{'xBIOMASS', 'y001654'},[0 (-0.47*3)/CNratios(i)]);
    ecModel = setParam(ecModel,'obj','xBIOMASS',1);
    sol = solveLP(ecModel);
    
    ecModel = setParam(ecModel,'eq','xBIOMASS',-sol.f*0.99);
    ecModel = setParam(ecModel, 'obj', 'prot_pool_exchange', 1);
    sol   = solveLP(ecModel);
    printFluxes(ecModel,sol.x)

    ecModel = setParam(ecModel, 'lb', 'prot_pool_exchange', sol.x(strcmpi(ecModel.rxns, 'prot_pool_exchange')) * 1.01);
    ecModel = setParam(ecModel,'obj','EXC_OUT_m1640',1);
    sol = solveLP(ecModel);

    growthRates(i) = sol.x(growthPos);
    lipidProduction(i) = -sol.f;
    %printFluxes(ecModel,sol.x)
end

% Create a figure for Growth Rate
figure;
plot(CNratios, growthRates, 'o-', 'LineWidth', 2);
xlabel('consumed CN Ratio');
ylabel('Growth Rate');
title('Growth Rate vs CN Ratio');
grid on;
xlim([0, max(CNratios)]);

% Create a separate figure for Lipid Production
figure;
plot(CNratios, lipidProduction, 'o-', 'LineWidth', 2);
xlabel('consumed CN Ratio');
ylabel('Lipid Production');
title('Lipid Production vs CN Ratio');
grid on;
xlim([0, max(CNratios)]);
