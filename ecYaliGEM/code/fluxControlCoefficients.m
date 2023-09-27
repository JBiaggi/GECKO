adapterLocation = fullfile(findGECKOroot,'ecYaliGEM','ecYaliGEMAdapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();

ecModel = loadEcModel('ecYaliGEM_FSEOF_pooled.yml');

ecModel = setParam(ecModel,'obj','xBIOMASS',1);

sol = solveLP(ecModel,1);

ecModel = setParam(ecModel,'lb','xBIOMASS',-sol.f*0.99);

% Fix carbon source and protein pool uptakes
ecModel = setParam(ecModel,'eq','y001808',sol.x(strcmp(ecModel.rxns,'y001808')));

ecModel = setParam(ecModel,'obj','EXC_OUT_m1640',1);

ecModel_OG = ecModel;

sol = solveLP(ecModel,1);
ecModel = setParam(ecModel,'eq','prot_pool_exchange',sol.x(strcmp(ecModel.rxns,'prot_pool_exchange')));

lipOG = -sol.f;

Qa = 1.01;

FCCs = zeros(numel(ecModel.ec.kcat),1);

% TO-DO: fix protein pool

progressbar('Flux control coefficients calculation')
for i = 1:numel(ecModel.ec.kcat)
    ecModel.ec.kcat(i) = ecModel.ec.kcat(i)*Qa;
    ecModel = applyKcatConstraints(ecModel);
    sol = solveLP(ecModel,1);
    Qlip = -sol.f;
    %FCCs(i) = 1000*(Qlip-lipOG)/lipOG;
    kcat_i = ecModel_OG.ec.kcat(i);
    FCCs(i) = (Qlip - lipOG)*kcat_i/(lipOG*(kcat_i * Qa - kcat_i));
    ecModel = ecModel_OG;
    progressbar(i/numel(ecModel.ec.kcat))
end
progressbar(1)

% Create output variable 
result.FCCs = FCCs;
result.kcat = ecModel_OG.ec.kcat;
result.rxns = ecModel_OG.ec.rxns;

% Organize output variable
[result.FCCs, sortOrder] = sort(result.FCCs,'descend');
result.kcat = result.kcat(sortOrder);
result.rxns = result.rxns(sortOrder);

for i = 1:length(result.rxns)
    result.rxnNames(i,1) = ecModel.rxnNames(strcmp(ecModel.rxns,result.rxns(i)));
end

