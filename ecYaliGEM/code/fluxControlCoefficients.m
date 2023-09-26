adapterLocation = fullfile(findGECKOroot,'ecYaliGEM','ecYaliGEMAdapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();

ecModel = loadEcModel('ecYaliGEM_FSEOF_pooled.yml');

ecModel = setParam(ecModel,'obj','xBIOMASS',1);

sol = solveLP(ecModel,1);

ecModel = setParam(ecModel,'lb','xBIOMASS',-sol.f*0.99);

ecModel = setParam(ecModel,'eq','y001808',sol.x(strcmp(ecModel.rxns,'y001808')));

ecModel = setParam(ecModel,'obj','EXC_OUT_m1640',1);

ecModel_OG = ecModel;

sol = solveLP(ecModel,1);

lipOG = -sol.f;

Qa = 1.001;

FCCs = zeros(numel(ecModel.ec.kcat),1);

% TO-DO: fix protein pool

progressbar('Flux control coefficients calculation')
for i = 1:numel(ecModel.ec.kcat)
    ecModel.ec.kcat(i) = ecModel.ec.kcat(i)*Qa;
    ecModel = applyKcatConstraints(ecModel);
    sol = solveLP(ecModel,1);
    Qlip = -sol.f;
    FCCs(i) = 1000*(Qlip-lipOG)/lipOG;
    ecModel.ec.kcat = ecModel_OG.ec.kcat;
    progressbar(i/numel(ecModel.ec.kcat))
end
progressbar(1)